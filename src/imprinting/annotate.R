#!/usr/bin/env Rscript
#
# Annotate DMRs based on other genomic features (including
# allele-specific H3K27me3 at promoters).
#
# run("imprinting", "annotate")
dependOn("imprinting", c("define_regions", "topol"))

# Make a big annotation table with all different types of information about DMRs and DEGs combined:
curGenes <- geneAnnotationExpr$gene
curRegsDt <- regsDt[imprintedRegionCandidates[ranges_type=="dmrseq",regionId],]
curRegs <- lib$dtToGr(curRegsDt,"seqnames")
curRegsM <- lib$dtToGr(curRegsDt[imprint_type=="M",],"seqnames")
curRegsP <- lib$dtToGr(curRegsDt[imprint_type=="P",],"seqnames")
dtAll <- geneAnnotationExpr
dtAll <- merge(dtAll, dtExprComplete, by="gene", all=F)[!is.na(chrom),]

dmrInPromo <- overlapsAny(lib$dtToGr(geneAnnotationExpr[,.(chrom, start=tss-promoWin, end=tss+promoWin)]), lib$dtToGr(curRegsDt,"seqnames"), ignore.strand=T)
dmrInGbody <- overlapsAny(lib$dtToGr(geneAnnotationExpr[,.(chrom, start, end)]), lib$dtToGr(curRegsDt,"seqnames"), ignore.strand=T) & !dmrInPromo
dmrUpstream <- overlapsAny(lib$dtToGr(geneAnnotationExpr[,.(chrom, start=ifelse(strand>0, tss, tss-100000), end=ifelse(strand>0, tss+100000, tss))]), lib$dtToGr(curRegsDt,"seqnames"), ignore.strand=T) & !dmrInPromo & !dmrInGbody
dmrDownstream <- overlapsAny(lib$dtToGr(geneAnnotationExpr[,.(chrom, start=ifelse(strand>0, start-100000, end), end=ifelse(strand>0, start, end+100000))]), lib$dtToGr(curRegsDt,"seqnames"), ignore.strand=T) & !dmrInPromo & !dmrInGbody

dtAll[, has_dmr_in_promo := dmrInPromo]
dtAll[, has_dmr_in_gbody := dmrInGbody]
dtAll[, has_dmr_100kb_us := dmrUpstream]
dtAll[, has_dmr_100kb_ds := dmrDownstream]
dtAll[, has_dmr_100kb_usds := dmrUpstream | dmrDownstream]

dtAll[, is_confirmed:=gene%in%geneLists[list_type=="is_confirmed",gene]]
dtAll[, is_unconfirmed:=gene%in%geneLists[list_type=="is_unconfirmed",gene]]
dtAll[,is_bix_up := is_bix & log2fc>0]
dtAll[,is_bix_down := is_bix & log2fc<0]
dtAll[,is_bsx_up := is_bsx & log2fc>0]
dtAll[,is_bsx_down := is_bsx & log2fc<0]


### add allele-specific H3K27me3 peak annotation ###

# load data obtained from GEO:
loadLibraries(c("R.utils","liftOver"))
ch <- resultsDir("mm9toMm10.chain.gz")
if(!file.exists(gsub(".gz","",ch))) {
	msg("download mm9 to mm10 chain")
	download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz", destfile=ch)
	gunzip(ch)
}
ch <- import.chain(gsub(".gz","",ch))
f <- resultsDir("GSE76687_ICM_k27me3_maternal_broadpeak.bed.gz")
if(!file.exists(f)) download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE76687&format=file&file=GSE76687%5FICM%5Fk27me3%5Fmaternal%5Fbroadpeak%2Ebed%2Egz", destfile=f)	
h3k27me3mat <- fread(f)
h3k27me3mat <- lib$dtToGr(h3k27me3mat,"V1","V2","V3")
h3k27me3mat <- unlist(liftOver(h3k27me3mat, ch))
f <- resultsDir("GSE76687_ICM_k27me3_paternal_broadpeak.bed.gz")
if(!file.exists(f)) download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE76687&format=file&file=GSE76687%5FICM%5Fk27me3%5Fpaternal%5Fbroadpeak%2Ebed%2Egz", destfile=f)	
h3k27me3pat <- fread(f)
h3k27me3pat <- lib$dtToGr(h3k27me3pat,"V1","V2","V3")
h3k27me3pat <- unlist(liftOver(h3k27me3pat, ch))
h3k27me3Regs <- c(h3k27me3mat,h3k27me3pat)
# - annotate all genes with a H3K27me3 peak within 5kb of their TSS:
curDist <- 5000
curGeneWindows <- lib$dtToGr(transcriptAnnotationExpr[,.(chrom, start=tss-curDist, end=tss+curDist)])
dtAll[, h3k27me3_imprint:=gene%in%transcriptAnnotationExpr[overlapsAny(curGeneWindows, h3k27me3Regs, ignore.strand=T), unique(gene)]]

# load H3K27me3 data from gametes:
f <- resultsDir("GSE76687_h3k27me3_sperm.broadpeak.gz")
if(!file.exists(f)) download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2041066&format=file&file=GSM2041066%5FSperm%5Fk27me3%5Fbroadpeak%2Ebed%2Egz", destfile=f)	
h3k27me3sp <- fread(f)
h3k27me3sp <- lib$dtToGr(h3k27me3sp,"V1","V2","V3")
h3k27me3sp <- unlist(liftOver(h3k27me3sp, ch))

f <- resultsDir("GSE76687_h3k27me3_oocyte.broadpeak.gz")
if(!file.exists(f)) download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76687/suppl/GSE76687_MII_oocyte_k27me3_broadpeak.bed.gz", destfile=f)	
h3k27me3oo <- fread(f)
h3k27me3oo <- lib$dtToGr(h3k27me3oo,"V1","V2","V3")
h3k27me3oo <- unlist(liftOver(h3k27me3oo, ch))

f <- resultsDir("GSE76687_h3k27me3_icm_any.broadpeak.gz")
if(!file.exists(f)) download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE76687&format=file&file=GSE76687%5FICM%5Fk27me3%5Fbroadpeak%2Ebed%2Egz", destfile=f)	
h3k27me3icmany <- fread(f)
h3k27me3icmany <- lib$dtToGr(h3k27me3icmany,"V1","V2","V3")
h3k27me3icmany <- unlist(liftOver(h3k27me3icmany, ch))

# annotate genes overlapping with H3K27me3 peaks:
dtAll[, h3k27me3_icm_any:=gene%in%transcriptAnnotationExpr[overlapsAny(curGeneWindows, h3k27me3icmany, ignore.strand=T), unique(gene)]]
dtAll[, h3k27me3_icm_mat:=gene%in%transcriptAnnotationExpr[overlapsAny(curGeneWindows, h3k27me3mat, ignore.strand=T), unique(gene)]]
dtAll[, h3k27me3_icm_pat:=gene%in%transcriptAnnotationExpr[overlapsAny(curGeneWindows, h3k27me3pat, ignore.strand=T), unique(gene)]]
dtAll[, h3k27me3_sperm:=gene%in%transcriptAnnotationExpr[overlapsAny(curGeneWindows, h3k27me3sp, ignore.strand=T), unique(gene)]]
dtAll[, h3k27me3_oocyte:=gene%in%transcriptAnnotationExpr[overlapsAny(curGeneWindows, h3k27me3oo, ignore.strand=T), unique(gene)]]

#m <- lib$dtToDf(dtAll[,.(gene, is_bix, is_bix_up, is_bix_down, is_bsx, is_bsx_up, is_bsx_down, h3k27me3_icm_mat, h3k27me3_icm_pat, h3k27me3_imprint, h3k27me3_sperm, h3k27me3_oocyte)]) * 1
	
# define sperm- or oocyte-specific H3K27me3 as those seen in one but not the other:
dtAll[, h3k27me3_sperm_specific:=h3k27me3_sperm&!h3k27me3_oocyte]
dtAll[, h3k27me3_oocyte_specific:=!h3k27me3_sperm&h3k27me3_oocyte]
# likewise for maternal/paternal in the ICM:
dtAll[, h3k27me3_icm_mat_specfic:=h3k27me3_icm_mat&!h3k27me3_icm_pat]
dtAll[, h3k27me3_icm_pat_specfic:=h3k27me3_icm_pat&!h3k27me3_icm_mat]

# define H3K27me3 imprinting status: 2 or -2 denotes sperm/oocyte/ICM-maternal/ICM-paternal specific, 1 denotes unspecific, 0 denotes no H3K27me3 present at promoter:
dtAll[, h3k27me3_sperm_status:=ifelse(h3k27me3_sperm,ifelse(h3k27me3_sperm_specific, 2, 1),0)]
dtAll[, h3k27me3_oocyte_status:=ifelse(h3k27me3_oocyte,ifelse(h3k27me3_oocyte_specific, -2, -1),0)]
dtAll[, h3k27me3_gamete_status:=ifelse(h3k27me3_oocyte|h3k27me3_sperm,ifelse(h3k27me3_oocyte_specific, -2, ifelse(h3k27me3_sperm_specific, 2, 1)),0)]
dtAll[, h3k27me3_icm_status:=ifelse(h3k27me3_icm_mat_specfic,-2,ifelse(h3k27me3_icm_pat_specfic,2,ifelse(h3k27me3_icm_any|h3k27me3_icm_pat|h3k27me3_icm_mat, 1, 0)))]

# generate a heatmap summarizing H3K27me3 patterns:

m <- lib$dtToDf(dtAll[,.(gene, is_bix_down=is_bix_down*3, is_bix_up=is_bix_up*3, is_bsx_down=is_bsx_down*3, is_bsx_up=is_bsx_up*3,h3k27me3_gamete_status, h3k27me3_icm_status)]) * 1
m <- m[rowSums(m[,grepl("is_b[si]x",colnames(m))]!=0)>0,]
a <- lib$dtToDf(dtAll[rownames(m),.(gene, dir=sign(log2fc), gene_list=ifelse(is_confirmed==T, "is_confirmed", ifelse(is_nbix==T, "is_nbix", "is_nbsx")))])
aCols <- list(dir=brewer.pal(3,"RdBu"), gene_list=colorPalettes$genesetLabel[unique(a[,2])]) 
aCols$gene_list[["is_nbsx"]] <- alpha(aCols$gene_list[["is_nbsx"]], 0.25)
m <- m[hclust(dist(m))$order,]
m <- m[order(a[rownames(m),"dir"], a[rownames(m),"gene_list"]),]
	colnames(m) <- gsub("h3k27me3_(.+)_status","\\1",colnames(m))
pheatmap::pheatmap(m[,5:6], gaps_row=cumsum(table(a[rownames(m),1])), annotation_colors=aCols, annotation_row=a, cluster_cols=F, cluster_rows=F, main="All BsX", breaks=-3.5:3.5, file=resultsDir("hm_h3k27me3_oo_sp_status_filt.pdf"), treeheight_col=0, treeheight_row=0, cellheight=9, cellwidth=10, col=c("black","red","yellow","white","yellow","blue","black"))
pheatmap::pheatmap(m[,5:6], gaps_row=cumsum(table(a[rownames(m),1])), show_rownames=F, annotation_colors=aCols, annotation_row=a, cluster_cols=F, cluster_rows=F, main="All BsX", breaks=-3.5:3.5, file=resultsDir("hm_h3k27me3_oo_sp_status_filt_nolbl.pdf"), treeheight_col=0, treeheight_row=0, cellheight=1, cellwidth=10, col=c("black","red","yellow","white","yellow","blue","black"))
fwrite(cbind(a[rownames(m),], m[,5:6]), file=resultsDir("source_data_figure_4a.csv"))

# genenerate pie charts:
ref <- "is_robust_expr"
pData <- melt(melt(dtAll, measure.vars=c(ref, "is_bix","is_bsx", "is_confirmed", "is_unconfirmed", "is_nbix", "is_nbsx", "published_imprint", "hc_repositories_imprint", "is_equivalent_0.1", "h3k27me3_confirmed_inoue_2017", "h3k27me3_candidate_inoue_2017", "h3k27me3_imprint"), variable.name="gene_selection", value.name="v1")[v1==T,], measure.vars=c("h3k27me3_gamete_status", "h3k27me3_icm_status"), variable.name="dev_stage")[, .(gene_selection=ifelse(gene_selection==ref,as.character(gene_selection),paste_(gene_selection, sign(log2fc))), dev_stage, status=ifelse(value==-2, "M", ifelse(value==2, "P", ifelse(abs(value)==1, "both", "-"))))]	
pData <- pData[, .(.N),by=.(dev_stage, gene_selection, status)]
pData <- pData[, .(N2=N, N=sum(N), perc=N/sum(N), status), by=.(dev_stage, gene_selection)]
selCats <- c("is_nbix", "is_nbsx", "is_confirmed")
selCats <- c(ref,c(paste_(selCats,1),paste_(selCats,-1)))
p <- ggplot(pData[gene_selection%in%selCats,], aes(x="", y=perc, fill=status))+ geom_bar(width=1, stat="identity") + coord_polar("y", start=0) + theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"), axis.text.x=element_blank(), legend.title=element_blank()) + facet_wrap(gene_selection~dev_stage+N, dir="v",nrow=2)
p <- p + scale_fill_manual(values=c("-"="#EEEEEE", "both"="yellow", "M"="red", "P"="blue"))
gg(p, "pie_h3k27_devstage", 44, 6, type="pdf")

fwrite(pData[gene_selection%in%grep(ref,selCats,value=T),.(gene_selection, dev_stage, status, perc, N)], file=resultsDir("source_data_figure_4b.csv"))
fwrite(pData[gene_selection%in%grep("_-1$",selCats,value=T),.(gene_selection, dev_stage, status, perc, N)], file=resultsDir("source_data_figure_4c.csv"))
fwrite(pData[gene_selection%in%grep("_1$",selCats,value=T),.(gene_selection, dev_stage, status, perc, N)], file=resultsDir("source_data_figure_4d.csv"))

pvals <- rblapply(setdiff(pData[,unique(gene_selection)],ref), function(geneSel) {
	rblapply(pData[,as.character(unique(dev_stage))], function(devStage) {
		rblapply(c("M","P"), function(imprintDir) {
			msgF("%s %s %s", geneSel, devStage, imprintDir)
			tab <- lib$dtToDf(dcast(unique(pData[gene_selection%in%c(ref,geneSel) & dev_stage==devStage, .(N2=sum(N2),N), by=.(gene_selection, dev_stage, status=ifelse(status==imprintDir,imprintDir,paste0("X",imprintDir)))]), status~gene_selection, fun.aggregate=sum, value.var="N2"))#[,c(2,1)]
			print(tab)
			data.table(fg=tab[1,1]/tab[2,1], bg=tab[1,2]/tab[2,2], pval=fisher.test(tab, alternative="greater")$p.value)
		}, "imprint_dir")
	}, "dev_stage")
}, "gene_selection")
pvals[,padj:=p.adjust(pval)]
pvals[,psig:=lib$pToSig(padj)]
pvals[grepl("bix",gene_selection),]
fwrite(pvals, file=resultsDir("d_h3k27me3_pvals_fisher.csv"))

pData[, stat:=as.factor(status)]
res <- rblapply(c("is_bix","is_bsx"), function(geneSel) {
	res1 <- rblapply(c("1","-1"), function(geneDir) {
		data.table(fixed_size="gene_dir", pval=chisq.test(apply(lib$dtToDf(dcast(pData[gene_selection==paste_(geneSel,geneDir),], status~dev_stage, value.var="N2")), 1, lib$naToZero))$p.value)
	}, "fixed_val")
	res2 <- rblapply(c("h3k27me3_gamete_status","h3k27me3_icm_status"), function(devStage) {
		data.table(fixed_size="dev_stage", pval=chisq.test(apply(lib$dtToDf(dcast(pData[grepl(geneSel,gene_selection) & dev_stage==devStage,], status~gene_selection, value.var="N2")), 1, lib$naToZero))$p.value)
	}, "fixed_val")
	rbind(res1,res2)
}, "gene_selection")
res[, psig:=lib$pToSig(pval)]
fwrite(res, file=resultsDir("d_h3k27me3_pvals.csv"))




### annotate genes by distance to DMR ###

# annotate genes with next DMR or DMRs within a certain window:
dtAll[, best_overlap:=ifelse(h3k27me3_imprint==T,"H3K27me3","")]
dtAll[, best_dist:=ifelse(h3k27me3_imprint==T,-2,Inf)]
dtAll[, best_overlap_M:=best_overlap]
dtAll[, best_dist_M:=best_dist]
dtAll[, best_overlap_P:=best_overlap]
dtAll[, best_dist_P:=best_dist]
d <- distanceToNearest(lib$dtToGr(dtAll), curRegs)
dtAll[queryHits(d), c("dmr_region_id", "dmr_meth_dif", "dmr_padj", "dmr_group"):=curRegsDt[subjectHits(d), .(rid, meth_diff, padj, grp)]]
dtAll[queryHits(d), dmr_dist:=d@elementMetadata$distance]
for(curDist in distThreshs) {
	cname <- paste_("meth_gene",round(curDist/1000),"kb")
	lbl <- sprintf("Meth:Gene%dkb",curDist/1000)
	
	curGeneWindows <- lib$dtToGr(dtAll[,.(chrom, start=start-curDist, end=end+curDist)])
	
	dtAll[, (cname):=overlapsAny(curGeneWindows, curRegs)]
	dtAll[, (paste_(cname,"M")):=overlapsAny(curGeneWindows, curRegsM)]
	dtAll[, (paste_(cname,"P")):=overlapsAny(curGeneWindows, curRegsP)]
	
	dtAll[(curDist==0 | h3k27me3_imprint==F) & !grepl("meth",best_overlap,ignore.case=T) & get(cname)==T, c("best_dist","best_overlap"):=.(curDist,ifelse(best_overlap=="",lbl,paste0(best_overlap,"; ",lbl)))]
	dtAll[(curDist==0 | h3k27me3_imprint==F) & !grepl("meth",best_overlap_M,ignore.case=T) & get(paste_(cname,"M"))==T, c("best_dist_M","best_overlap_M"):=.(curDist,ifelse(best_overlap_M=="",lbl,paste0(best_overlap_M,"; ",lbl)))]
	dtAll[(curDist==0 | h3k27me3_imprint==F) & !grepl("meth",best_overlap_P,ignore.case=T) & get(paste_(cname,"P"))==T, c("best_dist_P","best_overlap_P"):=.(curDist,ifelse(best_overlap_P=="",lbl,paste0(best_overlap_P,"; ",lbl)))]
}

# annotate with membership in different gene lists:
for(lt in geneLists[,unique(list_type)]) {
	dtAll[, (lt):=gene%in%geneLists[list_type==lt,gene]]
}
dtAll[topol==T & !grepl("Meth",best_overlap), best_overlap:=ifelse(h3k27me3_imprint==T,"H3K27me3, same TAD","Same TAD")]

# genes are associated with at least one DMR if the DMR is within 250kb of the gene and/or if there's a DMR in the same TAD:
dtAll[, meth_imprint:=topol | (is.finite(dmr_dist) & dmr_dist<=250000)]

# plot pie charts for genes in different categories and segments indicating the distance to the next DMR or H3K27me3 imprint:
pData <- melt(dtAll, measure.vars=names(genesetLabels))[value==T,.(N=.N,genes=paste(sort(gene),collapse="; ")), by=.(variable,best_dist,best_overlap)]
pData[, list_type:=factor(variable, levels=names(genesetLabels), labels=genesetLabels)]
pData[, list_size:=lib$dt2namedVec(geneLists[,.N,by=.(ID=list_type)])[as.character(variable)]]
pData[,best_overlap:=factor(best_overlap,levels=pData[order(-best_dist),unique(best_overlap)])]
pData <- pData[order(best_overlap),]
pData[, perc:=N/sum(N), by=list_type]
pData[best_overlap=="", best_overlap:="> 250kb, outside TAD"]
plotCols <- c(
	"H3K27me3"="#cb181d",
	"H3K27me3; Meth:Gene0kb"="#6a51a3",
	"H3K27me3, same TAD"="#c994c7",
	"Meth:Gene0kb"="#253494",
	"Meth:Gene100kb"="#2c7fb8",
	"Meth:Gene250kb"="#41b6c4",
	"Same TAD"="#a1dab4",
	"> 250kb, outside TAD"="#efefec"
)
pData[, best_overlap:=factor(best_overlap, levels=names(plotCols))]
p <- ggplot(pData, aes(x="", y=perc, fill=best_overlap))+ geom_bar(width=1, stat="identity") + coord_polar("y", start=0) + theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"), axis.text.x=element_blank(), legend.title=element_blank()) + facet_wrap(~sprintf("%s\n(%s)\nn = %d",list_type,variable,list_size), ncol=8)
p <- p + scale_fill_manual(values=plotCols)
gg(p, "pie", 16, 8, type="pdf")

fwrite(pData[list_type%in%c("allexpressed", "published_imprint", "hc_repositories_imprint", "is_confirmed", "is_unconfirmed", "is_equivalent_0.1", "is_nbix", "is_nbsx"), .(perc, best_overlap, list_type, variable, list_size)][order(list_type),], file=resultsDir("source_data_figure_3b.csv"))
