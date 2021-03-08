#!/usr/bin/env Rscript
#
# Load input data (raw WGBS and processed RNA-seq), gene lists, etc.
#
# run("imprinting", "load")
dependOn("imprinting", c("init"))


### Load processed RNA-seq data: ###

dtExprComplete <- fread(dataDir("rnaseq_300320.csv"))

# redefine BiX / BsX labels based on earlier labelling (BiX = BiX_HS, BsX = BiX_HS | BiX_LS):
dtExprComplete$is_bix <- dtExprComplete$is_bix_hs
dtExprComplete$is_bsx <- dtExprComplete$is_bix_hs | dtExprComplete$is_bix_ls
dtExprComplete$is_nbix <- dtExprComplete$is_nbix_hs
dtExprComplete$is_nbsx <- dtExprComplete$is_nbix_hs | dtExprComplete$is_nbix_ls
dtExprComplete$is_nbsx_not_nbix <- !dtExprComplete$is_nbix & dtExprComplete$is_nbsx
dtExprComplete$is_bsx_not_bix <- !dtExprComplete$is_bix & dtExprComplete$is_bsx
dtExprComplete$is_bix_ls <- NULL
dtExprComplete$is_bix_hs <- NULL
dtExprComplete$is_nbix_ls <- NULL
dtExprComplete$is_nbix_hs <- NULL


### Define gene lists: ###

selCols <- c(
	"h3k27me3_confirmed_inoue_2017",
	"h3k27me3_candidate_inoue_2017",
	"published_imprint",
	"is_equivalent_0.1",
	"hc_repositories_imprint",
	"is_bix",
	"is_nbix",
	"is_bsx",
	"is_nbsx",
	"is_nbsx_not_nbix"
)

geneLists <- melt(dtExprComplete, measure.vars=selCols, id.vars="gene")[value>0,]
geneLists <- split(geneLists[,as.character(gene)], f=geneLists$variable)
geneLists$is_confirmed <- intersect(geneLists$published_imprint,geneLists$is_bsx)
geneLists$is_unconfirmed <- setdiff(geneLists$published_imprint,geneLists$is_bsx)
geneLists$allexpressed <- dtExprComplete$gene
geneLists <- sapply(sapply(geneLists, unique), sort)

allExprGenes <- geneLists$allexpressed

geneLists <- as.data.table(melt(geneLists))[,.(list_type=L1,gene=value)]
geneLists[,gene:=as.character(gene)]
geneLists[,.(.N, length(unique(gene))),by=list_type]


### Get annotation for all genes: ###

transcriptAnnotation <- getTranscriptAnnotation(-1)
# Collapse transcripts to genes:
geneAnnotation <- transcriptAnnotation[,.(id=unique(gene), strand=mean(strand), chrom=unique(chrom), start=min(start), end=max(end), tss=ifelse(mean(strand)>0, max(end), min(start))),by=gene]
setkey(geneAnnotation, id)
## n.b. in this annotation 1=reverse, -1=forward strand

# Refine the gene annotation to only include those that are robustly detectable (and detected) in our dataset.
# We define this here as those genes with at least a certain number of reads covering an unambiguous part of 
# the gene's exons, i.e. a part that is not overlapping with any other gene.
rnaseqBams <- grep("STAR",list.files(lib$baseResultsDir("rnaseq"), ".bam$", full.names=T), invert=T, value=T)
rnaseqNames <- gsub(".*/(.*).bam","\\1",rnaseqBams)
excludedSamples <- grepl("^cast.+bla_[46]$",rnaseqNames)
rnaseqBams <- rnaseqBams[!excludedSamples]
rnaseqNames <- rnaseqNames[!excludedSamples]
# N.B. we use simpleCache to save results of lengthy analysis steps (https://cran.r-project.org/web/packages/simpleCache/index.html)
simpleCache(paste_("robust_transcripts",robustTransReads,robustTransSamples,digest::digest(rnaseqNames)), {
	allExprNonPred <- allExprGenes
		
	# get all exons of expressed/detectable genes:
	ex <- getExonAnnotation()[gene%in%allExprNonPred,]
	ext <- 100 # also exclude "near-overlaps"
	extEx <- lib$dtToGr(ex[,.(chrom, start=start-ext, end=end+ext)])

	# create disjoint unique exon fragments:
	exFrags <- disjoin(extEx, ignore.strand=T)
	o <- findOverlaps(extEx, exFrags, ignore.strand=T)
	o <- data.table(ex_frag=subjectHits(o), ex[queryHits(o), .(gene)]) #[ ,.(.N,length(unique(gene)),ex_frag,gene), by=ex]

	# then remove all fragments that overlap with more than one gene exon:
	o <- o[!ex_frag%in%o[,.(N=length(unique(gene))),by=ex_frag][N>1, ex_frag],]	
	disjointExons <- unique(data.table(gene=o$gene, lib$grToDt(exFrags[o$ex_frag]), strand=0))

	# determine coverage in unambiguous exons per sample:
	loadLibrary("Rsubread")
	exCov <- Rsubread::featureCounts(rnaseqBams, useMetaFeatures=T, annot.ext=as.data.frame(disjointExons[,.(GeneId=gene, Chr=chrom, Start=start, End=end, Strand=strand)]), allowMultiOverlap=T, strandSpecific=0, nthreads=16)
	colnames(exCov$counts) <- rnaseqNames
	exCovX <- melt(as.data.table(exCov$counts, keep.rownames="gene"), id.vars="gene")

	fwrite(exCov$counts, file=resultsDir("d_disjointexons_counts.csv"))

	# get list of all genes matching criteria:
	robustExprTranscripts <- exCovX[value>=robustTransReads, .N, by=gene][N>=robustTransSamples, unique(gene)]

	robustExprTranscripts
}, noload=F, assignToVar="robustExprTranscriptsComputed", cacheDir=analysisCacheDir)
	
# take list of robust transcripts from table instead (previously added to input table):
robustExprTranscripts <- dtExprComplete[is_robust_expr==T, gene]

# the two lists (from table and calculated de-novo should be identical (because the table has been filled in with values 
# from a previous run of the same code). Perform a sanity check:
if(!assertthat::are_equal(sort(robustExprTranscripts),sort(robustExprTranscriptsComputed))) {
	stop("mismatch between annotated robust genes in input table and those computationally determined")
}
rm(allExprGenes)

# restrict gene lists to robust ones:
transcriptAnnotationExpr <- transcriptAnnotation[gene%in%robustExprTranscripts,]
geneAnnotationExpr <- geneAnnotation[gene%in%robustExprTranscripts,]
geneLists <- geneLists[gene%in%robustExprTranscripts,]





### Load DNA methylation data: ###

# Prepare cache for CpG methylation data (unfiltered):
simpleCache("cpg_meth", {
	dMeth <- rblapply(dA$sample_name, function(s) {
		msg(s)
		suppressWarnings(readMethFile(dA[s,paste("zcat",meth_file)]))
	}, "sample_name")
	dMeth
}, assignToVar="cpgMeth", noload=TRUE, cacheDir=analysisCacheDir)

# Summary of WGBS coverage etc.:
simpleCache("cpg_meth_summary", {
	simpleCache("cpg_meth", assignToVar="cpgMeth", cacheDir=analysisCacheDir)
	pData <- cpgMeth[sample_name%in%dA[use==1,sample_name],]
	pData <- pData[, .(.N, N_thresh=sum(readCount>=minReadsThresh), perc_meth=sum(hitCount)/sum(readCount)), by=sample_name]
	pData
}, assignToVar="cpgSummary", noload=F, cacheDir=analysisCacheDir)
cpgSummary[, lbl:=dA[sample_name, sample_label_short]]
p <- ggplot(suppressWarnings(melt(cpgSummary, id.vars="lbl", measure.vars=c("N", "perc_meth"))), aes(x=factor(lbl, levels=cpgSummary[order(perc_meth),lbl]), y=value)) + geom_bar(stat="identity") + defTheme(flipX=T) + facet_wrap(~variable, ncol=1, scales="free_y")
gg(p, "cpg_stats", 4, 3.5, type="pdf", addBase=T, expand0=T)
fwrite(suppressWarnings(melt(cpgSummary, id.vars="lbl", measure.vars=c("N", "perc_meth"))), file=resultsDir("source_data_extended_data_figure_2d.csv"))

# Aggregate CpG methylation data in commonly covered positions and convert to BSseq object:
simpleCache(paste_("all_raw", minSamplesThresh, minReadsThresh), {
	simpleCache("cpg_meth", assignToVar="dMeth", cacheDir=analysisCacheDir)

	dPos <- dMeth[readCount>=minReadsThresh,.N,by=.(chr,start)][order(chr,start),]
	dPos <- dPos[N>=minSamplesThresh,] # only keep positions that have been counted more than minSamplesThresh times with at least minReadsThresh reads
	dPos[, pos_id:=1:.N]

	dMeth <- merge(dMeth, dPos, by=c("chr","start"))[,.(pos_id, hits=hitCount, total=readCount, sample_name)]

	M <- dcast(dMeth, pos_id~sample_name, value.var="hits")
	Cov <- dcast(dMeth, pos_id~sample_name, value.var="total")

	M <- lib$dtToDf(M)[dPos$pos_id, dA$sample_name]
	Cov <- lib$dtToDf(Cov)[dPos$pos_id, dA$sample_name]
	M[is.na(M)] <- 0
	Cov[is.na(Cov)] <- 0

	bs <- BSseq(chr = dPos$chr, pos = dPos$start, M = as.matrix(M), Cov = as.matrix(Cov),  sampleNames = dA$sample_name)
	bs
}, assignToVar="bs", noload=TRUE, cacheDir=analysisCacheDir)
