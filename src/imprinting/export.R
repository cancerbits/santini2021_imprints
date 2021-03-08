#!/usr/bin/env Rscript
#
# Export complete annotation tables for DMRs and DEGs and
# generate plots comparing both.
#
# run("imprinting", "export")
dependOn("imprinting", c("annotate","expression","region_enrichments","external","motifs"))



### COMPLETE CANDIDATE REGION ANNOTATION ###

fc <- lib$dt2namedVec(dtExprComplete[,.(log2fc, ID=gene)])
	
msg("Generate complete novel imprint DMR annotation table:")
dtRegs <- regsDtCand[,.(region_id=rid,imprint_type,chrom=seqnames,start,end,width,n_cpgs=L,pval,padj,min_cov=minCov,total_cov=totalCov,beta_embryo_andro=beta_Embryo_andro,beta_embryo_partheno=beta_Embryo_partheno)]
dtRegs$meth_diff <- (rowMeans(dWide[dtRegs$region_id,dA[sample_group=="Embryo_andro",sample_name]]) - rowMeans(dWide[dtRegs$region_id,dA[sample_group=="Embryo_partheno",sample_name]])) * 100
gr <- lib$dtToGr(dtRegs)

# overlaps with known imprints
msg("* annotating overlaps with known imprints..")
o <- findOverlaps(gr, grFeil)
o <- data.table(id=queryHits(o), dtFeil[subjectHits(o),.(k)])[, .(.N, paste(k,collapse=",")), by=id]
dtRegs[, overlaps_known_num:=0L]
dtRegs[o$id, c("overlaps_known_num", "overlaps_known") := o[,.(N,V2)]]

# next gene + dist
msg("* annotating genes..")
o <- as.data.table(distanceToNearest(gr, lib$dtToGr(geneAnnotation)))
dtRegs[o$queryHits, nearest_gene_dist:=o$distance]
dtRegs[o$queryHits, nearest_gene:=geneAnnotation[o$subjectHits, gene]]

# next expressed gene + dist + expr:
msg("* annotating expressed genes..")
o <- as.data.table(distanceToNearest(gr, lib$dtToGr(geneAnnotationExpr)))
dtRegs[o$queryHits, nearest_expressed_gene_dist:=o$distance]
dtRegs[o$queryHits, nearest_expressed_gene:=geneAnnotationExpr[o$subjectHits, gene]]
dtRegs[o$queryHits, nearest_expressed_log2fc:=fc[geneAnnotationExpr[o$subjectHits, gene]]]

# next published imprinted gene:
msg("* annotating published imprinted genes..")
i <- geneAnnotationExpr$gene %in% geneLists[list_type=="published_imprint", gene]
gtmp <- geneAnnotationExpr[i,]
o <- as.data.table(distanceToNearest(gr, lib$dtToGr(gtmp)))
dtRegs[o$queryHits, nearest_published_imprinted_gene_dist:=o$distance]
dtRegs[o$queryHits, nearest_published_imprinted_gene:=gtmp[o$subjectHits, gene]]
dtRegs[o$queryHits, nearest_published_imprinted_gene_log2fc:=fc[gtmp[o$subjectHits, gene]]]

# next our imprinted genes:
for(sel in c("bsx","nbsx","bix","nbix")) {
	msg("* annotating ",sel," genes..")
	i <- geneAnnotationExpr$gene %in% geneLists[list_type==paste_("is", sel), gene]
	gtmp <- geneAnnotationExpr[i,]
	o <- as.data.table(distanceToNearest(gr, lib$dtToGr(gtmp)))
	dtRegs[o$queryHits, paste_("nearest",sel,"dist"):=o$distance]
	dtRegs[o$queryHits, paste_("nearest",sel):=gtmp[o$subjectHits, gene]]
	dtRegs[o$queryHits, paste_("nearest",sel,"log2fc"):=fc[gtmp[o$subjectHits, gene]]]
}

# motif overlaps
msg("* annotating motif overlaps..")
fimo <- rblapply(selMotifDbs, function(motifDb) {
	motifFile <- motifDbs[[motifDb]]$meme
	resultsFile <- resultsDir("fimo_out/",inputs$fimo,"_", motifDb, "_t0.001.txt.gz")
	fimo <- data.table(readFIMOResult(resultsFile, motifFile=motifFile, motifPThresh=motifPThresh))	
	fimo[, regionId:=inputs$dt[as.numeric(substring(fimo$seq,4)), regionId]]
	fimo
}, "motifDb")[motifName%in%motifsSelFocus, ]
for(n in motifsSelFocus) {
	univIds <- fimo[motifName==n,regionId]	
	o <- overlapsAny(gr,  lib$dtToGr(dtUniv[regionId%in%univIds,]))
	dtRegs[, paste0("motif_",n) := o]
}

# .. and add cluster IDs to DMR annotation:
dtRegs[, wang_cluster_id:=wangClustIds[region_id]] 

msg("* writing output..")
fwrite(dtRegs, file=resultsDir("tab_dmrs.csv"))




##### COMPLETE CANDIDATE GENE ANNOTATION #####

setkey(dtRegs, region_id)

msg("Generate complete novel imprinted gene annotation table:")
dtGenes <- dtAll[,.(gene_symbol=gene, chrom, start, end, strand, nearest_dmr=dmr_region_id, nearest_dmr_dist=dmr_dist, nearest_dmr_padj=dmr_padj, nearest_dmr_meth_diff=dmr_meth_dif, our_meth_imprint=meth_imprint)]

# allele-specific expression
msg("* annotating expression..")
dtGenes <- cbind(dtGenes, dtAll[, setdiff(colnames(dtExprComplete),c("gene","topol")), with=F]) 

# imprinting KOs
msg("* annotating imprinting KOs..")
dtGenes <- cbind(dtGenes, dtAll[, setdiff(grep("chen|inoue",colnames(dtAll),value=T), colnames(dtGenes)), with=F])

# H3K27me3-imprint reannotation
msg("* annotating H3K27me3-imprinted genes..")
dtGenes <- cbind(dtGenes, dtAll[, setdiff(grep("h3k27",colnames(dtAll), value=T), colnames(dtGenes)), with=F])


dtGenes <- cbind(dtGenes, dtAll[, setdiff(names(genesetLabels), colnames(dtGenes)), with=F]) 

msg("* writing output..")
fwrite(dtGenes, file=resultsDir("tab_genes.csv"))
fwrite(dtGenes[!is.na(padj_chen_ctrl) | !is.na(padj_inoue_ctrl),], file=resultsDir("tab_genes_filtfig5.csv"))







### Generate scatterplots contrasting allelic bias in gene expression with allelic bias in DMR methylation: ###

# Get gene sets and define colors:
selGeneCats <- c("is_nbsx", "is_nbix", "is_confirmed", "is_unconfirmed") #"is_equivalent_0.1", 
plotCols <- colorPalettes$genesetLabel[selGeneCats]
plotCols["-"] <- "grey"
plotShapes <- structure(rep(1,length(plotCols)), names=names(plotCols))
plotShapes[c("is_nbix", "is_confirmed")] <- 16
plotSizes <- structure(rep(1.5,length(plotCols)), names=names(plotCols))
plotSizes["-"] <- 0.1

# Make plots for genes vs. nearest DMR:
lfc <- log2(1.5)

pData <- dtGenes[abs(nearest_dmr_meth_diff)>=thresh_meth_diff & abs(log2fc)>=lfc,]
pData[,d:=nearest_dmr_dist]
pData[,gene_category:="-"]
for(lt in selGeneCats) {
	pData[gene_symbol%in%geneLists[list_type==lt, gene], gene_category:=lt] 
}
pData[, is_bg:= !gene_category%in%selGeneCats]

p1 <- ggplot(pData, aes(x=log2fc, y=nearest_dmr_meth_diff)) 
p1 <- p1 + geom_vline(xintercept=lfc, linetype="dashed") + geom_vline(xintercept=-lfc, linetype="dashed") + geom_hline(yintercept=thresh_meth_diff, linetype="dashed") + geom_hline(yintercept=-thresh_meth_diff, linetype="dashed")
p1 <- p1 + geom_point(aes(shape=gene_category, color=gene_category, size=gene_category, alpha=ifelse(is_bg, 0.25, 1)), alpha=1)
p1 <- p1 + theme_bw() + theme(legend.position="top", legend.title=element_blank()) + scale_alpha_identity() + scale_size_manual(values=plotSizes) + scale_color_manual(values=plotCols) + scale_shape_manual(values=plotShapes)
p1 <- p1 + geom_text_repel(aes(label=gene_symbol), size=2.5, force=3, point.padding=0.5, segment.alpha=0.5, min.segment.length=0, data=pData[is_confirmed | (is_bix & (our_meth_imprint)),]) + xlab("Delta RNA") + ylab("Delta DNA meth")
p1 <- p1 + coord_cartesian(xlim=c(-6,9.05), ylim=c(-100,100))
gg(p1, paste_("scatter_meth_expr",thresh_meth_diff,"250kborTAD","gs"), 6, 4, type="pdf")

fwrite(pData[, .(x=log2fc, y=nearest_dmr_meth_diff, gene_category, gene_symbol, show_label=is_confirmed | (is_bix & (our_meth_imprint)))], file=resultsDir("source_data_figure_3a.csv"))