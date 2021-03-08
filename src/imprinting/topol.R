#!/usr/bin/env Rscript
#
# Load data on TADs.
#
# run("imprinting", "topol")
dependOn("imprinting", c("define_regions"))

threshTopolScore <- 0.5
selTopolData <- 2:11

# load TADs (data from Collombet et al. (2020)):
dtTopol <- fread(dataDir("external/collombet_2020_tads.csv"))
dtTopol[, chrom:=gsub("^(.+)_(.+)_(.+)$","\\1",domain)]
dtTopol[, start:=as.numeric(gsub("^(.+)_(.+)_(.+)$","\\2",domain))]
dtTopol[, end:=as.numeric(gsub("^(.+)_(.+)_(.+)$","\\3",domain))]
grTopol <- lib$dtToGr(dtTopol)

grGenes <- lib$dtToGr(geneAnnotation[dtExprComplete$gene,])

for(topolCluster in dtTopol[, rev(unique(cluster))]) {
	grCur <- grTopol[dtTopol[, cluster==topolCluster], ]
	o <- overlapsAny(grGenes, grCur)
	dtExprComplete[, paste0("topol_clust_c",topolCluster):=o] # annotate TAD cluster overlaps
}

# now, only work with TADs that overlap a DMR:
o <- overlapsAny(grTopol, lib$dtToGr(imprintedRegionCandidates[ranges_type=="dmrseq",],"chr"))

dtTopol <- dtTopol[o,]
grTopol <- grTopol[o,]

dtTopol[, is_64c_tad:=`64CSE_G2`>=threshTopolScore | `64CSE_G1`>=threshTopolScore]

# annotate genes:
for(topolDataCol in colnames(dtTopol)[selTopolData]) {
	grCur <- grTopol[dtTopol[, get(topolDataCol)>=threshTopolScore], ]
	o <- overlapsAny(grGenes, grCur)
	dtExprComplete[, paste_("topol",topolDataCol):=o] # since we previously filtered to TADs overlapping with DMRs, this now denotes genes that overlap a TAD that also overlaps a DMR
}

# 64-cell-stage TADs best match our developmental stage:
dtExprComplete[, topol:= topol_64CSE_G1 | topol_64CSE_G2]

grTopolSel <- grTopol[dtTopol$is_64c_tad,]
o1 <- findOverlaps(grTopolSel, lib$dtToGr(imprintedRegionCandidates[ranges_type=="dmrseq",],"chr"))
o2 <- findOverlaps(grTopolSel, lib$dtToGr(transcriptAnnotationExpr))
fwrite(unique(cbind(dtTopol[is_64c_tad==T,][queryHits(o2),.(tad_hit_i=queryHits(o2),domain)], transcriptAnnotationExpr[subjectHits(o2),.(gene)])), file=resultsDir("tad_genes.csv"))

topologicalData <- unique(merge(data.table(tad_id=queryHits(o1), region_id=imprintedRegionCandidates[ranges_type=="dmrseq",regionId][subjectHits(o1)]), data.table(tad_id=queryHits(o2), gene=transcriptAnnotationExpr$gene[subjectHits(o2)]), by="tad_id", allow.cartesian=T)[,.(gene,region_id)])

selCols <- unique(c(selCols, "topol"))
geneLists <- unique(rbind(geneLists, data.table(
	list_type="topol", gene=dtExprComplete[topol==TRUE, gene]
)))
