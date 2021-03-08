#!/usr/bin/env Rscript
#
# Integrate with published DNA methylation data.
#
# run("imprinting", "external")
dependOn("imprinting", c("annotate","expression"))


### Wang et al. (2014) ###

simpleCache("wang2014_overlap", {
	simpleCache("wang2014", {
		# Load data (mm10):
		extData <- c(
			Oocyte="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386019/suppl/GSM1386019_oocyte_mc_CG_plus.bed.gz", 
			Sperm="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386020/suppl/GSM1386020_sperm_mc_CG_plus.bed.gz",
			Embryo_2cell="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386021/suppl/GSM1386021_2cell_mc_CG_plus.bed.gz",
			Embryo_2cell_mat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386021/suppl/GSM1386021_2cell_mc_CG_maternal_plus.bed.gz",
			Embryo_2cell_pat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386021/suppl/GSM1386021_2cell_mc_CG_paternal_plus.bed.gz",
			Embryo_4cell="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386022/suppl/GSM1386022_4cell_mc_CG_plus.bed.gz",
			Embryo_4cell_mat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386022/suppl/GSM1386022_4cell_mc_CG_maternal_plus.bed.gz",
			Embryo_4cell_pat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386022/suppl/GSM1386022_4cell_mc_CG_paternal_plus.bed.gz",
			Embryo_ICM="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386023/suppl/GSM1386023_ICM_mc_CG_plus.bed.gz",
			Embryo_ICM_mat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386023/suppl/GSM1386023_ICM_mc_CG_maternal_plus.bed.gz",
			Embryo_ICM_pat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386023/suppl/GSM1386023_ICM_mc_CG_paternal_plus.bed.gz",
			Embryo_E3.5="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2577nnn/GSM2577161/suppl/GSM2577161_E35_mc_CG_plus.bed.gz",
			Embryo_E6.5="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386024/suppl/GSM1386024_E65_mc_CG_plus.bed.gz",
			Embryo_E6.5_mat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386024/suppl/GSM1386024_E65_mc_CG_maternal_plus.bed.gz",
			Embryo_E6.5_pat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386024/suppl/GSM1386024_E65_mc_CG_paternal_plus.bed.gz",
			Embryo_E7.5="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386025/suppl/GSM1386025_E75_mc_CG_plus.bed.gz",
			Embryo_E7.5_mat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386025/suppl/GSM1386025_E75_mc_CG_maternal_plus.bed.gz",
			Embryo_E7.5_pat="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386025/suppl/GSM1386025_E75_mc_CG_paternal_plus.bed.gz",
			PGC_male_E13.5="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386027/suppl/GSM1386027_E135M_mc_CG_plus.bed.gz",
			PGC_female_E13.5="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1386nnn/GSM1386026/suppl/GSM1386026_E135F_mc_CG_plus.bed.gz"
		)
		dtExt <- rblapply(names(extData), function(n) {
			f <- resultsDir("ext_",n,".bed.gz")
			if(!file.exists(f)) download.file(extData[n], destfile=f)
			suppressWarnings(fread(paste("zcat",f)))
		}, "dataset")
		setnames(dtExt, c("dataset","chr","start","end","perc","n_meth","n_unmeth"))
		dtExt <- dtExt[!grepl("_",chr),]
		dtExt
	}, recreate=F, assignToVar="dtExt", cacheDir=analysisCacheDir)

	# Determine overlaps with DMRs:
	tmpDt <- imprintedRegionCandidates[ranges_type%in%c("dmrseq","feil"),]
	tmpGr <- lib$dtToGr(tmpDt,"chr")
	o <- findOverlaps(tmpGr, lib$dtToGr(dtExt,"chr"), ignore.strand=T)
	o2 <- data.table(regionId=tmpDt[queryHits(o),regionId], dtExt[subjectHits(o), .(n_meth,n_unmeth,perc,dataset)])
	o2
}, assignToVar="o2", cacheDir=analysisCacheDir)

# Plot heatmaps:
annot <- regsDtCand[,.(regionId=rid, lbl)]
setkey(annot, regionId)
m <- o2[, .(meth=sum(n_meth)/sum(n_meth+n_unmeth),meth2=mean(perc)), by=.(regionId, dataset)]
m <- lib$dtToDf(dcast(m, regionId~dataset, value.var="meth"))*100
mx <- m[feilOrder,c("Embryo_ICM", "Embryo_ICM_mat", "Embryo_ICM_pat", "Oocyte", "Sperm")]
pheatmap::pheatmap(mx, cellheight=10, cellwidth=10, file=resultsDir("heatmap_external_gl.pdf"), na_col="grey", col=methCols, breaks=seq(0,100,length.out=length(methCols)+1), show_colnames=T, show_rownames=T, cluster_rows=F, cluster_cols=F, border_color="white", main=sprintf("Known imprints in Wang 2014"))
fwrite(mx, file=resultsDir("source_data_extended_data_figure_2e.csv"))

lib$pdfPlot("heatmap_external",7, 4)
# All DMRs, only Oocyte/Sperm:
mx <- m[candRegsOrder,c("Oocyte","Sperm")]
# Keep clustering information / dendrogram:
wangHc <- hclust(dist(mx))
wangClustNameMap <- c(
	Wang_C1 = "DMR-C2",
	Wang_C2 = "DMR-C4",
	Wang_C3 = "DMR-C1",
	Wang_C4 = "DMR-C3",
	Wang_C5 = "DMR-C5"
	
)
wangClustIds <- as.character(wangClustNameMap[paste0("Wang_C",cutree(wangHc, k=5))])
names(wangClustIds) <- names(cutree(wangHc, k=5))
colorPalettes$dmrCluster <- c("#fc8d62", "#e78ac3", "#66c2a5", "#8da0cb", "#a6d854")
names(colorPalettes$dmrCluster) <- as.character(sort(wangClustNameMap))
wangHm <- pheatmap::pheatmap( mx[candRegsOrder,], annotation_colors=list(clust=colorPalettes$dmrCluster), annotation_row=lib$dtToDf(data.table(ID=rownames(mx), clust=as.character(wangClustIds))), cellwidth=10, na_col="white", cutree_row=5, show_rownames=F, col=methCols, breaks=seq(0,100,length.out=length(methCols)+1), show_colnames=T, cluster_rows=wangHc, cluster_cols=F, border_color="white", main=sprintf("Wang 2014 Sp/Oo"))
dev.off()

# Write table with cluster IDs:
fwrite(as.data.table(wangClustIds, keep.rownames="region_id"), file=resultsDir("wang_cluster_ids.csv"))


# Function to perform clustering only within pre-defined groups (strata).
# To be used to order regions with similar DNA methylation patterns next to each other but within
# the context of previously defined clusters.
stratifiedHclust <- function(x, strata, clustMethod="complete", distMethod="euclidean", verbose=F) {
	o <- c()	
	for(s in unique(strata)) {		
		i <- names(which(strata==s))
		if(verbose) message(sprintf("stratum '%s', %d members", s, length(i)))
		m <- x[i,]		
		hc <- hclust(dist(m, method=distMethod), method=clustMethod)
		o <- c(o, i[hc$order])
	}	
	o
}

lib$pdfPlot("heatmap_external_ext",7, 4)
mxExt <- cbind(mx, dWide[rownames(mx),rev(dA[grepl("Embryo",sample_group),][order(sample_group),sample_name])]*100)
pheatmap::pheatmap(mxExt, labels_row=annot[rownames(mxExt),lbl], gaps_col=2, annotation_colors=list(clust=colorPalettes$dmrCluster), annotation_row=lib$dtToDf(data.table(ID=rownames(mx), clust=as.character(wangClustIds))), cellwidth=10, na_col="white", cutree_row=5, col=methCols, breaks=seq(0,100,length.out=length(methCols)+1), show_colnames=T, cluster_rows=wangHc, cluster_cols=F, border_color="white", main=sprintf("Known imprints and embryo methylation in Wang 2014 Sp/Oo"))
pData <- mxExt[wangHc$order,]
colnames(pData)[-(1:2)] <- dA[colnames(pData)[-(1:2)],][, sample_label_short]
pData <- cbind(cluster=as.character(wangClustIds[rownames(pData)]), pData)
fwrite(pData, file=resultsDir("source_data_figure_2c.csv"))
mxExt3 <- cbind(mx, dWide[rownames(mx),rev(dA[grepl("Embryo",sample_group),][order(sample_group),sample_name])]*100, dWide[rownames(mx),rev(dA[!grepl("Embryo",sample_group),][order(sample_group),sample_name])]*100)
mxExt3 <- mxExt3[stratifiedHclust(mxExt3, strata=wangClustIds),]
tmp <- dA[colnames(mxExt3),sample_label_short]
pheatmap::pheatmap(mxExt3, labels_col=ifelse(is.na(tmp),colnames(mxExt3),tmp), labels_row=annot[rownames(mxExt3),lbl], gaps_row=cumsum(table(factor(wangClustIds[rownames(mxExt3)],levels=unique(wangClustIds[rownames(mxExt3)])))), gaps_col=c(2, 8, 9), annotation_colors=list(clust=colorPalettes$dmrCluster), annotation_row=lib$dtToDf(data.table(ID=rownames(mxExt3), clust=as.character(wangClustIds[rownames(mxExt3)]))), cellwidth=10, na_col="white", cutree_row=5, col=methCols, breaks=seq(0,100,length.out=length(methCols)+1), show_colnames=T, cluster_rows=F, cluster_cols=F, border_color="white", main=sprintf("Known imprints and all methylation in Wang 2014 Sp/Oo, reclustered within"))
pData <- mxExt3
colnames(pData)[-(1:2)] <- dA[colnames(pData)[-(1:2)],][, sample_label_short]
pData <- cbind(cluster=as.character(wangClustIds[rownames(pData)]), pData)
fwrite(pData, file=resultsDir("source_data_extended_data_figure_2k.csv"))
dev.off()

fwrite(as.data.table(wangClustIds, keep.rownames="region_id"), file=resultsDir("wang_newcluster_ids.csv"))

colorPalettes$regionLabelWang <- c(colorPalettes$regionLabel, colorPalettes$dmrCluster) 



### plot distribution of DMR clusters with respect to genes and CGIs ###

loadLibrary("annotatr")
	
simpleCache("annotatr_wang", {
	# Annotation categories (CGIs, enhancers, promoters, gene bodies):
	annotTypes <- c("mm10_cpg_islands") #"mm10_enhancers_fantom"
	regAnnots <- c(
		build_annotations(genome = genomeBuild, annotations = annotTypes),
		lib$dtToGr(unique(transcriptAnnotation[,.(chrom,start=tss-promoWin,end=tss+promoWin,id,type="promoter")]), metaCols=c("id","type")),
		lib$dtToGr(unique(geneAnnotation[,.(chrom,start,end,id,type="gene")]), metaCols=c("id","type"))
	)

	regsX <- rbind(
		imprintedRegionCandidates[ranges_type!="dmrseq",], 
		imprintedRegionCandidates[ranges_type=="dmrseq",.(chr, start, end, width, regionId, imprintType="wang", ranges_type="wang", grp=wangClustIds[regionId])],	
		randomTrials[,.(chr, start, end, width=end-start, regionId, imprintType="rnd", ranges_type="rnd", grp=i)]
	, fill=T)

	regsX[, grp:=gsub("_[MP]","",grp)]
	grpSizes <- lib$dt2namedVec(regsX[,.N,by=.(ID=grp)])

	# Intersect the tested regions with the annotations:
	tmpAnnotated <- rblapply(names(grpSizes), function(curGrp) {
		as.data.table(data.frame(annotate_regions(
			regions = lib$dtToGr(regsX[grp==curGrp,],"chr", metaCols=c("regionId")),
			annotations = regAnnots,
			ignore.strand = TRUE,
			quiet = FALSE
		)))
	},"grp")

	# Select interesting categories for simplicity:
	selCats <- c("cpg_islands","genes_introns","genes_intergenic","genes_promoters","mm10_basicgenes","mm10_enhancers_fantom")
	annotTypesShort <- gsub(paste0(genomeBuild,"_"),"",annotTypes)
	tmpAnnotated[,annot.type.short:= gsub(paste0(genomeBuild,"_"),"",annot.type)]
	tmpAnnotated[annot.type.short=="promoter" & !grepl("rnd",grp),length(unique(regionId)),by=grp]
	tmpAnnotated[, ridUniq:=paste_(regionId,start,end)] # multiple gene promoters have the same name and "regionId" --> make them unique here
	
	list(tmpAnnotated=tmpAnnotated, grpSizes=grpSizes)
}, recreate=F, assignToVar="tmpAnnotated", cacheDir=analysisCacheDir)
grpSizes <- tmpAnnotated$grpSizes
tmpAnnotated <- tmpAnnotated$tmpAnnotated

# plot results:
pData <- tmpAnnotated[,.(N=length(unique(ridUniq)),perc=length(unique(ridUniq))/grpSizes[as.character(grp)]),by=.(grp,annot.type.short)]
pData <- rbind(pData[,.(grp,annot=annot.type.short,perc,N,pos=T)], pData[,.(grp,annot=annot.type.short,N=NA,perc=1-perc,pos=F)])
pData[, metaGroup:=regionLabels[gsub("_[MP]$","",gsub("rnd_\\d+","rnd",grp))]]	
pData[is.na(metaGroup), metaGroup:=grp]
pData <- pData[pos==T,]

p <- ggplot(pData[pos==T & !grepl("rnd",grp),], aes(x=annot, y=perc*100)) + geom_violin(data=pData[pos==T & grepl("rnd",grp),], color=colorPalettes$regionLabel["Random regions"]) + defTheme(topLegend=T, noLegendTitle=T, flipX=F) + geom_point(aes(shape=metaGroup, col=metaGroup)) + xlab(NULL) + ylab("Regions overlapping genomic feature (%)") + coord_flip()
pLoc <- p + scale_color_manual(values=colorPalettes$regionLabelWang) + scale_shape_manual(values=sapply(pData[,unique(metaGroup)], function(x) ifelse(grepl("DMR-C",x),1,16)))
pdf(resultsDir("annot_distr_all_wang.pdf"), 6, 2.2, useDingbats=F)
print(pLoc)
dev.off()	

fwrite(pData[,.(grp,annot,perc,metaGroup)], file=resultsDir("source_data_figure_2d.csv"))
