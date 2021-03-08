#!/usr/bin/env Rscript
#
# Define imprinting DMRs.
#
# run("imprinting", "define_regions")
dependOn("imprinting", c("gold_standard","load"))

doAllPlots <- F
candRegsOrder <- NULL

# identify candidate differentially methylated regions (DMRs) using dmrseq:
simpleCache(paste_("dmrseq", "embryo_3v", minSamplesThresh, minReadsThresh, analysisVer), {
	simpleCache(paste_("all_raw", minSamplesThresh, minReadsThresh), assignToVar="bs", cacheDir=analysisCacheDir)

	pData(bs)[,groupVar] <- dA[,get(groupVar)]

	sampleI <- which(rownames(pData(bs)) %in% sampleSelections$embryo)
	bsFilt <- bs[, sampleI]
	lociI <- which(DelayedMatrixStats::rowSums2(getCoverage(bsFilt, type="Cov")==0) == 0)
	bsFilt <- bsFilt[lociI, ]

	pData(bsFilt)[,groupVar] <- factor(pData(bsFilt)[,groupVar], levels=c("Embryo_biparental","Embryo_andro","Embryo_partheno"))
		
	cl <- MulticoreParam(max(4,min(multicoreWorkers()/1.5, 16))) 

	msg("start dmrseq")
	regions <- dmrseq(bs=bsFilt, testCovariate=groupVar, BPPARAM=cl)		
	msg("finished dmrseq")
	
	tryCatch({stopCluster(cl)}, error=function(e) {return(TRUE)})
	regions
}, assignToVar="regions", reload=TRUE, recreate=recreateAll, cacheDir=analysisCacheDir)
regsDt <- as.data.table(data.frame(regions))
regsDt[, rid_num:=1:.N]
regsDt[, rid:=paste0("r", rid_num)]
regsDt[, ranges_type:="dmrseq"]
regsDt[, imprint_type:=ifelse(beta_Embryo_andro>0, "P", "M")]
regType <- "dmrseq"

# quantify read counts in the candidate regions across all samples:
simpleCache(paste_("reads_and_meth", regType, minSamplesThresh, minReadsThresh), {
	regSetX <- GRangesList(lib$dtToGr(regsDt,"seqnames"))
	names(regSetX) <- regType
	counts <- countAggregatedReads(dA, regSetX, name=paste_("rcountnmean",regType, minSamplesThresh, minReadsThresh), repeats=NULL, sizes=sapply(regSetX,length))
	counts
}, recreate=recreateAll, assignToVar="regReadCounts", reload=TRUE, cacheDir=analysisCacheDir)
regReadCounts[,regionId:=paste0("r",regionID)]

# add counts to DMR annotation table: 
regsDt <- merge(regsDt, regReadCounts[id %in% sampleSelections$embryo, .(minCov=min(readCount,na.rm=T), totalCov=sum(readCount,na.rm=T)), by=regionId], by.x="rid", by.y="regionId", all.x=T)
dWide <- lib$dtToDf(dcast(regReadCounts, regionId~id, value.var="methyl"))
regsDt[, meth_diff:=(rowMeans(dWide[rid,dA[sample_group=="Embryo_andro",sample_name]]) - rowMeans(dWide[rid,dA[sample_group=="Embryo_partheno",sample_name]])) * 100]

# apply independent pre-filtering criteria for regions to qualify as imprint DMRs:
regsDt[, opposing:=sign(beta_Embryo_andro * beta_Embryo_partheno)==-1] # check that DNA methylation difference is in diametrically opposing direction with respect to bi-parental control
regsDt[, prefilt:=opposing & totalCov>=minTotalCov & width>=100]
regsDt[prefilt==T, padj:=p.adjust(pval, "fdr")]
regsDt[, sig:=abs(beta_Embryo_andro)>=0.25 & abs(beta_Embryo_partheno)>=0.25 & abs(meth_diff)>=thresh_meth_diff & opposing==T & is.finite(padj) & padj<=0.1 ]	
regsDt[sig==T,.N]
regsDt[,grp:=paste_(ranges_type, imprint_type)]
regsDt[, dif:= beta_Embryo_partheno - beta_Embryo_andro]
regsDt[, dir:= sign(beta_Embryo_partheno)]

# all regions considered == background for enrichment analysis later on:
regionBg <- regsDt[prefilt==T,.(chr=seqnames, start, end, width, regionId=rid)]

# check for overlaps with known GL-DMRs:
o <- countOverlaps(grFeil,lib$dtToGr(regsDt,"seqnames"))	
regsDt[, nOverlapsFeil:=countOverlaps(lib$dtToGr(regsDt,"seqnames"),grFeil)]

# filter to imprint DMRs
regsDtCand <- regsDt[sig==T,]	
setkey(regsDtCand, rid)

simpleCache(paste_("feil", "regsDtExt"), {
	simpleCache(paste_("all_raw", minSamplesThresh, minReadsThresh), assignToVar="bs", cacheDir=analysisCacheDir)
	bsX <- unique(data.table(chrom=as.character(seqnames(bs)), start=start(bs), end=end(bs)))
	L <- countOverlaps(grFeil, lib$dtToGr(bsX))
	regsDtExt <- cbind(dtFeil, L=L, width=dtFeil$end-dtFeil$start)
	regsDtExt
}, recreate=F, assignToVar="regsDtExt", reload=TRUE, cacheDir=analysisCacheDir)

imprintedRegionCandidates <- regsDtCand[, .(chr=seqnames, start, end, width, regionId=rid, imprintType=imprint_type, ranges_type, nOverlapsFeil, grp)]

o <- findOverlaps(lib$dtToGr(regsDtCand,"seqnames"),grFeil,ignore.strand=T)
regsDtCand[queryHits(o), lbl:=dtFeil[subjectHits(o),name]]


# Aggregate allelic read coverage in bins around the center of all DMRs:
win <- 100000
binSize <- 100
nBins <- 100
simpleCache(paste_("dmrseq_agg",win, minSamplesThresh, minReadsThresh), {
	simpleCache("cpg_meth", assignToVar="cpgMeth", cacheDir=analysisCacheDir)

	cpgMethAvg <- rblapply(dA[sample_name%in%sampleSelections$embryo, unique(sample_group)], function(curGrp) {
		cpgMeth[sample_name%in%dA[sample_group==curGrp,sample_name],.(m=mean(hitCount/readCount)),by=.(chr,pos=start)]
	}, "grp")
			
	d <- rblapply(imprintedRegionCandidates$regionId, function(rid) {
		curReg <- imprintedRegionCandidates[regionId==rid,]
		center <- curReg[,round((end+start)/2)]
		
		curCpgs <- cpgMethAvg[chr==curReg[,chr] & pos>=center-win & pos<=center+win, ]
		curCpgs[, relPos:=pos-center]
		
		m <- dcast(curCpgs, relPos~grp, value.var="m", fun.aggregate=mean)[order(relPos),]
		m[,dif:=Embryo_andro-Embryo_partheno]
		
		m[,.(relPos,dif)]
	},"regionId")
	
	d
}, recreate=F, assignToVar="d", cacheDir=analysisCacheDir)
d[, relPosApprox:=round(relPos/binSize)]

rowAnnot <- lib$dtToDf(imprintedRegionCandidates[,.(regionId,imprintType)])
hmData <- lib$dtToDf(dcast(d[abs(relPos)<=binSize*nBins, ], regionId~relPosApprox, value.var="dif", fun.aggregate=mean))
hmData <- hmData[order(rowAnnot[rownames(hmData),"imprintType"],rowMeans(hmData,na.rm=T)),]

regsDtCand[, candRnk:=1:.N]
regsDtCand[is.na(lbl), lbl:=""]
setkey(regsDtCand, rid)
hmData <- hmData[order(rowAnnot[rownames(hmData),"imprintType"],regsDtCand[rownames(hmData),candRnk]),]

# impute missing values (for visualization only):
hmDataImp <- t(apply(hmData, 1, function(x) {
	i <- sample(which(is.na(x)))
	for(j in i) {
		lb <- max(0,j-2)
		ub <- min(j+2,ncol(hmData))
		x[j] <- mean(x[lb:ub],na.rm=T)
	}
	x
}))
dimnames(hmDataImp) <- dimnames(hmData)
		
# plot heatmap of aggregate coverage over bins in a window around all DMRs:
lib$pdfPlot("heatmap_dmrseq_agg", 6, 10)
pheatmap(hmDataImp, col=c("#550000",brewer.pal(9,"RdYlBu"),"#000055"), labels_row=regsDtCand[rownames(hmDataImp),lbl], gaps_row=cumsum(table(rowAnnot[rownames(hmData),"imprintType"])), na_col="white", show_colnames=F, show_rownames=T, cluster_rows=F, cluster_cols=F, annotation_row=rowAnnot, border_color="white", display_numbers=F, number_color="black", number_format="%.1f", annotation_colors=colorPalettes, main=sprintf("Candidate imprints, n = %d, bin size = %d bp, total width = %d kb, NAs imputed", nrow(hmData), binSize, nBins*binSize/1000))
dev.off()

# export source data file:
fwrite(
cbind(regsDtCand[rownames(hmDataImp),lbl], hmDataImp)
, file=resultsDir("source_data_figure_2b.csv"))

# remember the order of DMRs for other plots later on:
candRegsOrder <- rownames(hmDataImp)



# add reference/control regions to be used for comparison of DMRs later on:

# 1. promoter regions:
if(imprintedRegionCandidates[ranges_type=="ref",.N]==0) {

	promos <- unique(transcriptAnnotation[,.(chrom,start=tss-promoWin,end=tss+promoWin,gene)])
	
	# Only look at regions that "could have been detected" according to our pre-filtering criteria:
	o <- overlapsAny(lib$dtToGr(promos),lib$dtToGr(regionBg,"chr"))
	promos <- reduce(lib$dtToGr(promos[o,]))
	promos <- lib$grToDt(promos)	
	promos[, regionId:=paste_("promo_fragment",1:.N)]
	promos[, chr:=chrom]
	promos[, width:=end-start]

	imprintedRegionCandidates <- rbind(imprintedRegionCandidates,
		promos[, .(chr, start, end, width, regionId, imprintType="promoter", ranges_type="ref", nOverlapsFeil=0, grp=NA)]
	)
}

# 2. random region sets:
n <- imprintedRegionCandidates[ranges_type=="dmrseq",.N]
simpleCache(paste_("rndregs", "dmrseq", minSamplesThresh, minReadsThresh), {
	# Draw 1000 random samples of regions from the universe of all detectable regions:
	randomTrials <- rblapply(paste_("rnd",1:1000), function(i) {
		randomTrials <- regionBg[sample(nrow(regionBg),n),]
	}, "i")
}, recreate=F, assignToVar="randomTrials", cacheDir=analysisCacheDir)
if(sum(!overlapsAny(lib$dtToGr(imprintedRegionCandidates[ranges_type=="dmrseq",],"chr"), lib$dtToGr(regionBg,"chr")))>0) stop("sanity check failed: DMRs exist that are not in background")

# 3. known GL-DMRs (list from Feil lab):
if(imprintedRegionCandidates[ranges_type=="feil",.N]==0) {
	imprintedRegionCandidates <- rbind(imprintedRegionCandidates,
		dtFeil[, .(chr=chrom, start, end, width=end-start, regionId=name, imprintType=ifelse(paternal=="*","P","M"), ranges_type="feil", nOverlapsFeil=1, grp=NA)]
	)
}
imprintedRegionCandidates[is.na(grp), grp:=paste_(ranges_type, imprintType)]
forEach(imprintedRegionCandidates[,unique(ranges_type)], function(n) {
	o <- countOverlaps(grFeil,lib$dtToGr(imprintedRegionCandidates[ranges_type==n,],"chr"))
	msgF("%s: \t\t # overlaps candidates vs Feil = %d / %d, total # candidates = %d", n, sum(o>0), nrow(dtFeil), nrow(imprintedRegionCandidates[ranges_type==n,]))
})

# define the "universe" of all testable regions for enrichment analyses:
regUniv <- reduce(c(lib$dtToGr(regionBg,"chr"),lib$dtToGr(imprintedRegionCandidates,"chr")))
dtUniv <- lib$splitBigRegions(lib$grToDt(regUniv), 2500)
regUniv <- lib$dtToGr(dtUniv)
dtUniv[,k:=paste0("u",1:.N)]
for(rtype in imprintedRegionCandidates[,unique(ranges_type)]) {
	dtUniv[, paste(rtype):=overlapsAny(regUniv,lib$dtToGr(imprintedRegionCandidates[ranges_type==rtype,],"chr"))]
}
if(sum(!overlapsAny(lib$dtToGr(imprintedRegionCandidates,"chr"), regUniv))>0) stop("sanity check failed: DMRs exist that are not in background")
