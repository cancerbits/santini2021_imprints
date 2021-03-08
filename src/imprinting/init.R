#!/usr/bin/env Rscript
#
# Define constants and functions to be used throughout the analysis.
#
# run("imprinting", "init")

options(bitmapType='cairo')


# Set working directory and project name (assumes utility functions library is already loaded): 
lib$projectName <- "meth_imprint"
setwd(paste0(Sys.getenv("CODEBASE"),"/", lib$projectName))

# Load libraries:
loadLibraries(c("cowplot","TxDb.Mmusculus.UCSC.mm10.knownGene","annotatr", "matrixStats", "data.table", "ggplot2", "simpleCache", "pheatmap", "RColorBrewer", "viridis", "LOLA","ggsignif","ggrepel","dmrseq"))
tryCatch({library("Cairo")}, error=function(e) print(e))

analysisVer <- "santini_2021"
lib$setCurrentAnalysis(analysisVer)


# functions:

myHclust <- function(x, dist.meth="euclidean", clust.meth="complete") {
	dists <- lib$dist(x, method=dist.meth)
	dists[is.na(dists)] <- max(dists,na.rm=T)
	hclust(dists, method=clust.meth)
}
naTo <- function(x, to="") {
	 sapply(x, function(x) {
			if(is.null(x) | is.na(x) | length(x)==0) return(to)
			else return(x)
	})
}
runEnrichr <- function(geneLists, dbs=c("GO_Biological_Process_2018")) {
	enrichrRes <- rblapply(names(geneLists), function(curGrp) { #
		msg(curGrp)
		curGenes <- geneLists[[curGrp]]			

		res <- enrichr(curGenes, databases=dbs)
		rbindlist(res, "database", use.names=T, fill=T)
	}, "grp")

	enrichrRes[, n.hits:=as.numeric(gsub("^(\\d+)/(\\d+)$","\\1",Overlap))]
	enrichrRes[, n.total:=as.numeric(gsub("^(\\d+)/(\\d+)$","\\2",Overlap))]
	
	enrichrRes[, original.term:=Term]
	enrichrRes[, Term:=NULL]
	enrichrRes[, term:=sapply(original.term, function(x) {
		x <-  gsub(" (Mus musculus|Mouse)", " (mouse)", gsub(" (Homo sapiens|Human)", " (human)", gsub("_", " ", x), ignore.case=TRUE), ignore.case=TRUE)
		x <- gsub(" \\(NAFLD\\)", "", x)
		x <- gsub("^MP\\d+\\s+", "", x, perl=TRUE)
		x <- gsub("\\s+(hsa|WP)\\d+$", "", x, perl=TRUE)
		x <- gsub("\\s+\\(GO:.+\\)$", "", x, perl=TRUE)
		x <- gsub(" \\d+ ChIP-.+ ", " ", x, perl=TRUE)
		x <- gsub("positive", "pos.", x, perl=TRUE)
		x <- gsub("negative", "neg.", x, perl=TRUE)
		x <- gsub("regulation of", "reg. of", x, perl=TRUE)
		x <- gsub("involved in ", "in ", x, perl=TRUE)
		x <- gsub("ligase activity", "ligase", x, perl=TRUE)		
		x <- gsub("(GSE\\d+) sample \\d+", "\\1", x, perl=TRUE)
		x <- gsub("UMLS ", "", x)
		x <- gsub("\\s+", " ", x, perl=TRUE)
		x <- gsub("in DMSO-Rat-Primary rat ", "DMSO-Rat ", x)
		x <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "",x)
		x <- gsub(" (\\w+)*-\\w+-\\w+$","",x)
		x <- gsub(" .human.$","",x)
		x <- gsub(" GSE.+$","",x)
		x <- gsub(" (activity|kinase ARCHS4 coexpression)", "", gsub("(modification) process", "\\1", gsub("(tissue)? development", "dev.", gsub("(catabo|metabo)lic process","\\1ism",gsub("biosynthetic process","biosynthesis", gsub("([Ss]ignaling|regulatory) [Pp]athway", "pathway", gsub("\\s+$","",gsub(".human. (\\w+-?)*$","",x))))))))
		x <- gsub("GTEX-[^\\s]+ (.+ (fe)?male).+$","GTEX \\1",x)
		lib$capFirst(substr(x, 1, 64))
	})]
	enrichrRes[, uniq.term:=paste_(database,original.term)]

	enrichrRes
}
loadAnnot <- function() {
	f <- list.files(metaDir(), pattern="samples_", full.names=T)
	f <- f[sapply(f, file.size)>0] # drop empty files
	dA <- rblapply(f, fread, "annot_file")
	dA[, meth_file:=baseResultsDir("pipeline/", sample_name, "/bismark_", genomeBuild, "/extractor/", sample_name, ".aln.dedup.filt.CpG_report_filt.min.gz")]	
	setkey(dA, "sample_name")
	dA
}
readMethFile <- function(f, format=NULL) {
	if(is.null(format)) {
		format <- gsub(".+\\.(\\w+)$", "\\1", gsub(".gz$", "", tolower(f)), perl=TRUE)
	}
	dt <- fread(f)
	setnames(dt, c("chr", "start", "hitCount","readCount"))
	return(dt)
}
makeWide <- function(d, regSet) {	
	dWide <- as.data.frame(dcast(d, regionID~id, value.var="methyl"))
	rownames(dWide) <- dWide[,1]
	dWide <- as.matrix(dWide[,-1])
	i <- as.numeric(rownames(dWide))
	rownames(dWide) <- regSet[i, k]
	dWide <- unique(dWide)
	dWide
}
writeFASTA <- function(sequences, sequenceFile) {
	msg("\t* create input file: ", sequenceFile)
	tmp <- paste0(rep(">seq", length(sequences)*2), floor(1:(length(sequences)*2) / 2 +1))
	tmp[seq(2,length(tmp),by=2)] <- as.character(sequences)
	write.table(tmp, file=sequenceFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
	sequenceFile
}
readFIMOResult <- function(fimoResultFile, motifFile, nSeqs=NULL, motifPThresh=0.05) {
	cmd <- paste("grep 'MOTIF '",motifFile)
	msg(cmd)
	motifNames <- tryCatch(lib$dt2namedVec(fread(cmd, select=c(2,3), header=FALSE), "V2"), error=function(e) {
		tmp <- fread(cmd, select=c(2), header=FALSE)
		structure(tmp[,V2], names=tmp[,V2])
	})
	
	dt <- fread(fimoResultFile, select=c(1,3,8))
		
	if(is.null(nSeqs)) nSeqs <- dt[,length(unique(sequence_name))]
	dt[, motifName:=motifNames[motif_id]]

	dt[`p-value`<=motifPThresh, .(seq=sequence_name, motifName, motifId=motif_id)]
}
runFIMO <- function(sequences, outFile="fimo_out.txt", motifP=1e-4, motifFile=paste0(Sys.getenv("RESOURCES"),"/motifs/motif_databases_feb2017/JASPAR/JASPAR_CORE_2016_vertebrates.meme"), params="--no-qvalue --text --bgfile motif-file", fimoExec= toolsDir("meme/bin/fimo"), inputIsFile=length(sequences)==1) {

	dir.create(dirname(outFile), recursive=TRUE, showWarnings=FALSE)

	if(!inputIsFile) {
		sequenceFile <- writeFASTA(sequences, paste0(outFile, "_input.fasta"))
	}
	else {
		msg("\t* use existing input file: ", sequences)
		sequenceFile <- sequences
	}

	cmd <- paste(fimoExec,paste("--thresh",motifP), params,motifFile,sequenceFile," | gzip >",outFile)
	msg("\t* run external tool: $ ", cmd)
	system(cmd)

	if(!inputIsFile) {
		msg("\t* remove temporary input file: ", sequenceFile)
		file.remove(sequenceFile)
	}
	
	outFile
}
augmentLolaRes <- function(lolaResAll, qThresh=0.05, qvalCol="qvalMean", pvalCol="pvalMean", orCol="oddsRatioMean") {
	lolaResAll[, mLog10Q:=-10*log10(get(qvalCol))]
	lolaResAll[, mLog10P:=-10*log10(get(pvalCol))]

	lolaResAll[,log2odds:=log2(get(orCol))]

	lolaResAll[,qCap:=pmax(get(qvalCol), qThresh*0.001)]


	lolaResAll[,term:=antibody]
	lolaResAll[is.na(term),term:=gsub("_\\(.+\\)$","",gsub("GSM\\d+_","",gsub("Human_","",gsub("wgEncode.wg","",gsub(".(bed|narrowPeak)","",filename)))))]
	lolaResAll[collection=="sheffield_dnase", term:=paste0("DNase #",gsub(".bed","",filename), " (", sapply(strsplit(description,";"),function(x) paste(substr(x,1,3),collapse=";")), ")")]
	lolaResAll[,term:=gsub("^(EGFP|C)\\W","",gsub("_\\(.+\\)$","",toupper(term)))]

	
	# harmonize synonyms:
	lolaResAll[,term:=gsub("EP300", "P300", term)]
	lolaResAll[,term:=gsub("PU\\.?1", "SPI1", term)]
	lolaResAll[,term:=gsub("[−-]", "", term)]
	lolaResAll[,term:=gsub("POL(II|LL)", "POL2", term)]
	lolaResAll[,term:=gsub("POLIII", "POL3", term)]
	lolaResAll[antibody%in%c("GR","NR3C1"), antibody:="NR3C1"]
	lolaResAll[grepl("^E?P300",antibody,ignore.case=T), antibody:="EP300"]
	lolaResAll[antibody%in%c("GABP"), antibody:="GABPA"]
	lolaResAll[term=="P300",term:="EP300"]
	lolaResAll[term=="ERALPHA_A",term:="ESR1"]
	lolaResAll[term=="TCF7L2_C9B9",term:="TCF7L2"]

	# kick out Pol2:
	lolaResAll <- lolaResAll[!grepl("^POL",antibody,ignore.case=T),]

	lolaResAll
}
curateRegionDb <- function(regionDB) {
	capFirst <- function(str) paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))
	regionDB$regionAnno[,cellType:=capFirst(gsub("(cell|progenitor|precursor|phage|cyte|blast)s","\\1", tolower(cellType), perl=TRUE))]
	regionDB
}
getRefSeqFile <- function() {
	f <- resultsDir("refseq.txt.gz")
	if(!file.exists(f)) {
		download.file("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/ncbiRefSeq.txt.gz", destfile=f)
	}
	#	
	return(f)
}
getTranscriptAnnotation <- function(promo.size=500) {
	dt <- fread(getRefSeqFile(), select=c(2,4,13,5,6,3), col.names=c("id", "strand","gene", "start", "end", "chrom"))
	dt[grepl("Snrpn|Snurf",gene), gene:="Snurf/Snrpn"]	# since one isoform Snrpn is identitical with the annotation of Snurf, we find these indistinguishable and merge them into one gene annotation
	dt[, strand:=ifelse(strand!="+",1,-1)]
	dt[,tss:=start]
	dt[strand==1,tss:=end] # in this gene defintion, "-1" means forward, "1" means reverse strand!
	if(promo.size>=0) {
		dt[, start:=tss-promo.size]
		dt[, end:=tss+promo.size]
	}
	setkey(dt, id)
	dt
}
getExonAnnotation <- function() {
	dt <- fread(getRefSeqFile(), select=c(2,4,13,3, 10, 11), col.names=c("id","strand","gene", "chrom", "ex_starts", "ex_ends"))
	dt[grepl("Snrpn|Snurf",gene), gene:="Snurf/Snrpn"]		# since one isoform Snrpn is identitical with the annotation of Snurf, we find these indistinguishable and merge them into one gene annotation
	
	dt <- rbindlist(apply(dt, 1, function(x) {
		s <- strsplit(x[["ex_starts"]],",")[[1]]
		e <- strsplit(x[["ex_ends"]],",")[[1]]
		s <- as.numeric(s)
		e <- as.numeric(e)
		data.table(id=x[["id"]],gene=x[["gene"]], strand=x[["strand"]], chrom=x[["chrom"]], start=s, end=e)
	}))	
	dt[, strand:=ifelse(strand!="+",1,-1)]
	dt[strand==1,tss:=end] # in this gene defintion, "-1" means forward, "1" means reverse strand!	
	
	unique(dt)
}
# functions borrowed from Nathan Sheffield's genomics libraries:
countAggregatedReads <- function(dA, aggrRegions, name, repeats=NULL, sizes=sapply(aggrRegions,length), fileCol="meth_file") {
	cols=c("readCount", "methyl")
	funcs=c("sum", "mean")
	jCommand = buildJ(cols,funcs)

	dA[,meth_file_x:=get(fileCol)]
	dA[grepl(".gz",meth_file_x),methCallFileX:=paste("zcat",meth_file_x)]

	rbindlist(lapplyAlias(dA[,sample_name], sampleSummaryByRegion,
		regions=aggrRegions, excludeGR = repeats, regionsGRL.length = sizes, 
		dA, idColumn = "sample_name", fileColumn="meth_file_x", 
		cachePrepend=name,  cacheDir=analysisCacheDir, cacheSubDir=paste0("/rcount/", name), 
		jCommand=jCommand, byRegionGroup=F,
		readFunction=function(x) { tmp <- readMethFile(x); tmp[,methyl:=hitCount/readCount]; tmp }, recreate=FALSE, mc.preschedule=FALSE))
}
buildJ = function(cols, funcs) {
	r = paste("list(", paste(paste0(cols, "=", funcs, "(", cols, ")"), collapse=","), ")")
	return(r);
}
setLapplyAlias = function(cores=0) {
	if (cores < 1) {
		return(getOption("mc.cores"))
	}
	if(cores > 1) { #use multicore?
		if (requireNamespace("parallel", quietly = TRUE)) {
			options(mc.cores=cores);
   		} else {
			warning("You don't have package parallel installed. Setting cores to 1.")
			options(mc.cores=1); #reset cores option.
  		}
	} else {
		options(mc.cores=1); #reset cores option.
	}
}
lapplyAlias = function(..., mc.preschedule=TRUE) {
	if (is.null(getOption("mc.cores"))) { setLapplyAlias(1) }
	if(getOption("mc.cores") > 1) {
		return(parallel::mclapply(..., mc.preschedule=mc.preschedule))

	} else {
		return(lapply(...))
	}
}
nlist = function(...) {
 	fcall = match.call(expand.dots=FALSE)
	l = list(...);
	if(!is.null(names(list(...)))) { 
		names(l)[names(l) == ""] = fcall[[2]][names(l) == ""]
	} else {	
		names(l) = fcall[[2]];
	}
	return(l)
}
sampleSummaryByRegion = function(id, regions, excludeGR, annotationTable,
			 			idColumn, fileColumn, cachePrepend, cacheSubDir, 
						jCommand, readFunction, byRegionGroup= FALSE,
						regionsGRL.length=NULL, ...) {
	# Check 
 	if (! "GRanges" %in% class(regions) & !"GRangesList" %in% class(regions) & !"CompressedGRangesList" %in% class(regions) ) {
		 print(class(regions))
		stop("regions is not a GRanges");
	}

	# Make sure all inputs exist:
	inputs = list("id", "regions", "excludeGR", "annotationTable", "idColumn", "fileColumn", "cachePrepend", "cacheSubDir", "jCommand", "readFunction")
	isn = sapply(inputs, is.null)
	if (any(isn)) {
		message("null input")
		print(cbind(inputs, ex, isn))
	}

	if (! id %in% annotationTable[,get(idColumn)]) {
		stop("id ", id, " not in annotation table.")
	}

	# identify the sample:
	i = which (id == annotationTable[,get(idColumn)]); 
	cacheName = paste0(cachePrepend, id);
	infile = annotationTable[i, get(fileColumn)];
	message(i, ": ", cacheName);
	# Produce a cache; read in the file and summarize by regions
	simpleCache(cacheName, {
		inObject = readFunction(infile);
		inObject = BSAggregate(inObject, regions, excludeGR=excludeGR, regionsGRL.length=regionsGRL.length, jCommand=jCommand, byRegionGroup=byRegionGroup)
		if (is.null(inObject)) { 
			return(NULL)
		}
		inObject[, id:=id]
		inObject
	},
	cacheSubDir=cacheSubDir, 
	buildEnvir=nlist(infile, regions, jCommand, id, readFunction, excludeGR, regionsGRL.length, byRegionGroup), ...) # end simpleCache
	if (exists(cacheName)) {
		return(get(cacheName))
	} else {
		return(NULL)
	}
}
BSdtToGRanges = function(dtList) {
	gList = list();
	for (i in 1:length(dtList)) {
		#dt = dtList[[i]];
		setkey(dtList[[i]], chr, start)
		#convert the data into granges object
		gList[[i]] = GRanges(seqnames=dtList[[i]]$chr, ranges=IRanges(start=dtList[[i]]$start, end=dtList[[i]]$start), strand=rep("*", nrow(dtList[[i]])), hitCount=dtList[[i]]$hitCount, readCount=dtList[[i]]$readCount)
	}
	return(gList);
}
BSAggregate <- function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, sumCols=NULL, jCommand=NULL, byRegionGroup=FALSE, keep.na=FALSE) {

	if( "GRanges" %in% class(regionsGRL)) {
		regionsGRL = GRangesList(regionsGRL);
	} else if (! "GRangesList" %in% class(regionsGRL) && ! "CompressedGRangesList" %in% class(regionsGRL)) {
		stop("regionsGRL is not a GRanges or GRangesList object");
	}

	if(! is.null(excludeGR)) {
		BSDT = BSFilter(BSDT, minReads=0, excludeGR)
	}

	bsgr = BSdtToGRanges(list(BSDT));

	colModes = sapply(BSDT,mode);
	if (is.null(sumCols)) {
		sumCols = setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
		# Restrict to numeric columns.		
		sumCols = intersect(sumCols, names(colModes[which(colModes == "numeric")]))

	}
	regionsGR = unlist(regionsGRL)
	
	if(is.null(regionsGRL.length)) {
		if (length(regionsGRL) > 100) {
		message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF=FALSE)
		}
		regionsGRL.length = sapply(regionsGRL, length)
		message("Done counting regionsGRL lengths.");
	}

	# Build a table to keep track of which regions belong to which group
	region2group = data.table(
		regionID=1:length(regionsGR), 
		chr=as.vector(seqnames(regionsGR)), 
		start=as.vector(start(regionsGR)), 
		end=as.vector(end(regionsGR)),
		withinGroupID= unlist(sapply(regionsGRL.length, seq)),
		regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
	setkey(region2group, regionID)

	message("Finding overlaps...");
	fo = findOverlaps(bsgr[[1]], regionsGR)

	setkey(BSDT, chr, start)

	message("Setting regionIDs...");
	BSDT = BSDT[queryHits(fo),] #restrict the table to CpGs in any region.
	BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.

	if (is.null(jCommand)) {
		cols=c(sumCols, keepCols)
		funcs = c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
		jCommand = buildJ(cols, funcs)
	}
	message("jCommand: ", jCommand)
	
	# Define aggregation column. aggregate by region or by region group?
	if (byRegionGroup) {
		agCol = "regionGroupID";
	} else {
		agCol = "regionID"; # Default
	}

	# Build the by string
	if (is.null(splitFactor)) {
		byString = paste0("list(regionID)");
	} else {
		byString = paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
	}

	# Now actually do the aggregate:
	message("Combining...");
	bsCombined = BSDT[,eval(parse(text=jCommand)), by=eval(parse(text=byString))]
	setkey(bsCombined, regionID)

	# Define aggregation column. aggregate by region or by region group?
	if (byRegionGroup) {
		# must set allow=TRUE here in case there are multiple IDs (splitCol)
		bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
		if (is.null(splitFactor)) {
			byStringGroup = "list(regionGroupID)"
		} else {
			byStringGroup = paste0("list(", paste("regionGroupID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
		}
		
		bsCombined=bsCombined[,eval(parse(text=jCommand)), by=eval(parse(text=byStringGroup))]
		return(bsCombined);
	} else {
		e = region2group[bsCombined,]
		setkey(e, regionID);
		return(e);
	}
}
# from RnBeads:
rnb.beta2mval <- lib$rnb.beta2mval




# define constants:

genomeBuild <- "mm10"
genomeRef <- "BSgenome.Mmusculus.UCSC.mm10"
dA <- loadAnnot()[use==1,]

sampleSelections <- list(
	all=dA$sample_name,
	embryo=dA[grepl("embryo",sample_group,ignore.case=T),sample_name],
	ESC=dA[grepl("ESC",sample_group,ignore.case=T),sample_name]
)

groupVar <- "sample_group"

analysisCacheDir <- resultsDir("rcache")

colorPalettes <- list(
	"Methylation hub"=c("0"="white","1"="black")
	,sample_group=c(Kidney_somatic="#33a02c", Embryo_biparental="#6a3d9a", Embryo_partheno="#e31a1c", Embryo_andro="#1f78b4", ESC_andro_haploid="#a6cee3", ESC_partheno_haploid="#fb9a99", ESC_fused="#cab2d6", ESC_female_diploid="#fccde5", ESC_male_diploid="#8dd3c7")
	,sex = c("female"="#e41a1c", "male"="#377eb8")
	,flowcell = lib$getCategoryColors(dA[,unique(flowcell)], "Greys")
	,tissue = c("embryonic stem cell"="#4daf4a","embryo"="#ff7f00","kidney"="#a65628")
	,imprint_type = c("maternal"="#984ea3", paternal="#ffd92f")
)
methCols <- (colorRampPalette(brewer.pal(7, "YlGnBu"))(20)) # brewer.pal(5,"RdYlBu")

minSamplesThresh <- 2
minReadsThresh <- 3
promoWin <- 1000

dir.create(analysisCacheDir, showWarnings=F)

recreateAll <- FALSE

lib$setCurrentAnalysis(analysisVer)

distThreshs <- c(0, 100, 250) * 1000

regionLabels <- c(
	dmrseq_P = "Candidates (paternal)",
	dmrseq_M = "Candidates (maternal)",
	feil_P = "Known (paternal)",
	feil_M = "Known (maternal)",
	dmrseq = "Candidates",
	feil = "Known",
	ref_promoter = "Gene promoters",
	ref = "Gene promoters",
	rnd = "Random regions",
	all = "All"
)
colorPalettes$regionLabel <- c(
	"#daa420",
	"#daa620",
	"#00028b",
	"#00048b",
	"#daa520",
	"#00008b",
	"black",
	"black",
	"lightgrey",
	"black"
)
names(colorPalettes$regionLabel) <- regionLabels

colorPalettes$genesetLabel <- c(
	"published_imprint" = "#3182bd",
	"hc_repositories_imprint" = "#3182bd",
	"is_bsx" = "#756bb1",
	"is_bix" = "#756bb1",
	"is_bsx_not_bix" = "#756bb1",
	"is_confirmed" = "#00008b",
	"is_nbsx" = "#f6d109",
	"is_nbix" = "#f6d109",
	"is_nbsx_not_nbix" = "#f6d109",
	"allexpressed" = "darkgrey",
	"is_equivalent_0.1" = "#ff7f00",
	"is_unconfirmed" = "#cab2d6",
	"h3k27me3_candidate_inoue_2017"="#cb181d",
	"h3k27me3_confirmed_inoue_2017"="#cb181d"
)
names(genesetLabels) <- genesetLabels <- names(colorPalettes$genesetLabel)

robustTransReads <- 12
robustTransSamples <- 4

thresh_meth_diff <- 30
minTotalCov <- 100

doAllPlots <- F