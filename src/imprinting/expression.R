#!/usr/bin/env Rscript
#
# Generate plots related to the allele-specific RNA-seq analysis.
#
# run("imprinting", "expression")
dependOn("imprinting", c("load","annotate"))

# Get maternal- and paternal-specific read counts:
X <- data.matrix(as.data.frame(dtExprComplete[,grepl("^(forward|reverse)",colnames(dtExprComplete)),with=F]))
totReads <- colSums(rbind(X[,seq(1,ncol(X),by=2)],X[,seq(2,ncol(X),by=2)]))
X <- t(t(X)/totReads[sort(rep((1:(ncol(X)/2)), 2))])*1000000
mat <- grep("maternal", colnames(X))
pat <- grep("paternal", colnames(X))
		
# Plot heatmaps for the genes in all pre-defined gene lists:
exprHms <- list() # store heatmaps so that we can retrieve clustering order later on
gLists <- geneLists[, unique(list_type)]
pData <- rblapply(setdiff(gLists,"allexpressed"), function(gList) {

	gs <- geneLists[list_type==gList,as.character(gene)]
	gsKnown <- geneLists[list_type=="published_imprint",as.character(gene)]
	i <- dtExprComplete$gene%in%gs

	XX <- X[i,]
	rownames(XX) <- dtExprComplete$gene[i]

	hmD <- t(apply(log2(XX+1),1,scale))
	dimnames(hmD) <- dimnames(XX)
	
	# export source data files:
	if(gList=="is_bix") fwrite(hmD, file=resultsDir("source_data_figure_1a.csv"))

	exprHms[[gList]] <<- pheatmap::pheatmap(hmD, file=resultsDir("heatmap_expression_",gList,".pdf"), border_color="white", cellwidth=8, labels_row=ifelse(rownames(XX)%in%gsKnown,rownames(XX),""), cutree_row=2, cutree_col=2, scale="row", col=colorRampPalette(brewer.pal(5,"RdBu"))(21)) #

	XX <- t(apply(XX, 1, function(x) x/max(x)))
	data.table(
		gene=dtExprComplete$gene[i], 
		mat=rowMeans(XX[,mat]), 
		pat=rowMeans(XX[,pat]), 
		dif=rowMeans(XX[,pat])-rowMeans(XX[,mat])
	)

}, "gene_list")

# Plot parallel coordinates (maternal/paternal) for all gene lists:
pData2 <- melt(pData, measure.vars=c("mat","pat"))[!gene_list%in%c("allexpressed"),]

p <- ggplot(pData2, aes(x=variable, y=value, color=dif)) + geom_point() + defTheme() + facet_wrap(~gsub("_","\n",gsub("is_","",gene_list)), nrow=1)
p <- p + geom_line(aes(group=gene)) + scale_color_gradientn(colors=(brewer.pal(9,"RdBu")))
p <- p + geom_text(aes(label=gene), data=pData2[gene%in%intersect(geneLists[list_type=="equivalent",gene], geneLists[list_type=="unconfirmed",gene]) & variable=="pat" & gene_list=="unconfirmed",], color="black")
gg(p, "expr_par_coord", length(gLists)*1.3, 4, type="pdf")


# export source data files:
fwrite(dcast(pData2[gene_list%in%c("is_confirmed","is_unconfirmed","is_nbix"), ], gene_list+gene~variable, value.var="value"), file=resultsDir("source_data_figure_1b.csv"))
fwrite(dcast(pData2[gene_list%in%c("is_nbsx","is_nbsx_not_nbix"), ], gene_list+gene~variable, value.var="value"), file=resultsDir("source_data_extended_data_figure_1b.csv"))



# Define gene sets based on DMR/H3K27me3 imprinting:

geneLists <- geneLists[!list_type%in%c("h3k27me3_only", "dmr_only", "dmr_and_h3k27me3"), ]

geneLists <- rbind(geneLists, data.table(
	list_type="h3k27me3_only", gene=dtAll[!meth_imprint & h3k27me3_imprint,gene]
))
geneLists <- rbind(geneLists, data.table(
	list_type="dmr_only", gene=dtAll[meth_imprint & !h3k27me3_imprint,gene]
))
geneLists <- rbind(geneLists, data.table(
	list_type="dmr_and_h3k27me3", gene=dtAll[meth_imprint & h3k27me3_imprint,gene]
))

geneLists <- rbind(geneLists, data.table(
	list_type="is_nbsx_not_nbix", gene=dtAll[is_nbsx & !is_nbix, gene]
))				

gs <- rblapply(c("is_confirmed", "is_unconfirmed", "is_nbix", "is_nbsx_not_nbix"), function(gList) {
	dtExprComplete[gene%in%geneLists$gene[geneLists$list_type==gList], .(gene, log2fc)]
}, "list_type")[gene%in%robustExprTranscripts,]
setkey(gs, gene)
geneLists[,.N,by=list_type]
gs[,.N,by=list_type]
 



### Analyze gene lists in allele-specific scRNA-seq dataset (Borensztein et al.): ###
 
stageOrder <- c("Oo","2C","4C","8C","16C","32C","64C") # fixed order by dev. stage
	

# Load data:
f <- resultsDir("GSE80810_AllelicRatio.xls")
if(!file.exists(f)) {
	download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80810&format=file&file=GSE80810%5FAllelicRatio%2Exls%2Egz", destfile=paste0(f,".gz"))	
	gunzip(paste0(f,".gz"))
}
dExt <- cbind(as.data.table(readxl::read_xls(f, sheet=1))[,-2], as.data.table(readxl::read_xls(f, sheet=2))[,-(1:2)])
for(cname in colnames(dExt)[-1]) {
	dExt[, paste(cname):=as.numeric(get(cname))]
}

X <- data.matrix(lib$dtToDf(dExt))
XX <- as.data.frame(X)[gs[order(list_type,log2fc),gene],]
regex <- "^(KO_)?(Oo|\\d+C)_(CB|BC|B6|Cast)_\\d+.?$"
colAnnot <- data.frame(stage=gsub(regex,"\\2",colnames(XX)),genotype=gsub("^$","WT",gsub(regex,"\\1",colnames(XX))),cross=gsub(regex,"\\3",colnames(XX)))
rownames(colAnnot) <- colnames(XX)

gsColors <- colorPalettes$genesetLabel[gs[,unique(list_type)]] 
gsColors[4] <- alpha(gsColors[4],0.5)
names(gsColors) <- gs[,unique(list_type)]  

# Plot heatmap:
XX2 <- XX[, order(factor(colAnnot$stage,levels=stageOrder), factor(colAnnot$genotype,levels=c("WT","KO_")), colAnnot$cross)]
XX2 <- XX2[rowSums(is.finite(data.matrix(XX2)))>=3,]
logXX2 <- log2(XX2+0.5)

pheatmap::pheatmap(logXX2, 
	cluster_cols=F, cluster_rows=F, treeheight_row=0,
	breaks=seq(-1,1,length.out=6),
	main="Borensztein et al. (2017)", border_color="white", na_col="#FEFEFE",
	show_colnames=F, show_rownames=F,
	cellheight=1, cellwidth=1,  
	annotation_row=lib$dtToDf(gs[rownames(logXX2),.(gene,list_type,log2fc)]),
	annotation_col=colAnnot, 
	annotation_colors=list(
		log2fc=brewer.pal(5,"RdYlBu"),
		cross=structure(brewer.pal(4,"Paired"),names=c("B6","BC","Cast","CB")),
		stage=structure(plasma(length(stageOrder)),names=stageOrder),
		genotype=c("WT"="steelblue","KO_"="orange"),
		list_type=gsColors
	), 
	fontsize=6.5,
	file=resultsDir("heatmap_ext_borensztein_allelic_custom_ordered_compact.pdf"),
	cutree_row=2, scale="none", 
	gaps_col=cumsum(as.data.table(colAnnot[colnames(XX2),])[,.N,by=stage][,N]), 
	gaps_row=cumsum(gs[rownames(XX2),][,.N,by=list_type][,N]), 
	col=colorRampPalette(brewer.pal(9,"RdYlBu"))(5)
)

# export source data files:
pData <- lib$dtToDf(merge(
	lib$dtToDf(gs[rownames(logXX2),.(gene,list_type,log2fc)]), 
	rbind(t(colAnnot), t(apply(logXX2,1,function(x) sapply(x, function(y) sprintf("%.5f", y))))),
	by="row.names", all=T
))

# export source data file:
fwrite(pData[c(colnames(colAnnot), rownames(logXX2)), ], file=resultsDir("source_data_extended_data_figure_1j.csv"))









### Compare to maternal KO of Eed and Dnmt3l: ###

loadLibrary("DESeq2")
 
select <- "bsx"
minNtest <- 0 #20

pThresholds <- c(0.1, 0.01, 0.001)

# Get selected gene sets:
gs <- rblapply(c("is_confirmed", "is_nbix", "is_nbsx_not_nbix", "is_unconfirmed"), function(gList) {
	dtExprComplete[gene%in%geneLists$gene[geneLists$list_type==gList], .(gene, log2fc)]
}, "list_type")[gene%in%robustExprTranscripts,]
setkey(gs, gene)




#-- Load first maternal KO dataset: --#


f <- resultsDir("aay7246_Table_S1.xlxs")
if(!file.exists(f)) {
	download.file("https://advances.sciencemag.org/highwire/filestream/222532/field_highwire_adjunct_files/1/aay7246_Table_S1.xlsx", destfile=f)	
}
dExt <- as.data.table(readxl::read_xlsx(f, sheet=2))[-1,-c(2:5, 8, 11, 14, 17)]
colnames(dExt) <- c("V1", paste_(sort(rep(c("ctr1","ctr2","matko1","matko2"),2)), c("maternal","paternal")))
for(cname in colnames(dExt)[-1]) {
	dExt[, paste(cname):=as.numeric(get(cname))]
}
counts <- data.matrix(lib$dtToDf(dExt))
counts <- rbind(counts, "Snurf/Snrpn"=counts["Snurf",]) # we found Snurf/Snrpn indistinguishable (same genomic coordinates) and merged the genes in our annotations; use one of the count rows for comparison
# focus on comparable genes
counts <- counts[intersect(rownames(counts),dtExprComplete$gene),]

dAx <- data.frame(genotype=gsub("^(.+)\\d_.+$","\\1",colnames(counts)), allele=gsub("^.+\\d_(.).+$","\\1",colnames(counts)), rep=gsub("^(.+\\d)_.+$","\\1",colnames(counts)), stringsAsFactors=F)
rownames(dAx) <- colnames(counts)

# estimate size factors on merged MAT and PAT counts:
sampleCounts <- counts[,c(1,3,5,7)] + counts[,c(2,4,6,8)]
sampleDA <- data.frame(genotype=gsub("^(.+)\\d_.+$","\\1",colnames(sampleCounts)), stringsAsFactors=F)
ddsSample <- DESeqDataSetFromMatrix(countData = sampleCounts,
							  colData = sampleDA,
							  design = ~genotype)
ddsSample <- DESeq(ddsSample, betaPrior=F)

j <- rownames(counts)%in%gs$gene

# in control:
i <- dAx$genotype=="ctr"
ddsCTRL <- DESeqDataSetFromMatrix(countData = counts[j,i],
							  colData = dAx[i,],
							  design = ~rep+allele)
sizeFactors(ddsCTRL) <- sizeFactors(ddsSample)[ceiling(which(i)/2)]
ddsCTRL <- DESeq(ddsCTRL, betaPrior=F)
resCTRLallelic <- as.data.table(results(ddsCTRL, contrast=c("allele","m","p"), altHypothesis="greaterAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resCTRLallelic[padj<=pThresholds[1],]
resCTRLnonallelic <- as.data.table(results(ddsCTRL, contrast=c("allele","m","p"), altHypothesis="lessAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resCTRLnonallelic[padj<=pThresholds[1],]

# in mKO:
i <- dAx$genotype!="ctr"
ddsKO <- DESeqDataSetFromMatrix(countData = counts[j,i],
							  colData = dAx[i,],
							  design = ~rep+allele)
sizeFactors(ddsKO) <- sizeFactors(ddsSample)[ceiling(which(i)/2)]
ddsKO <- DESeq(ddsKO, betaPrior=F)
resKOallelic <- as.data.table(results(ddsKO, contrast=c("allele","m","p"), altHypothesis="greaterAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resKOallelic[padj<=pThresholds[1],]
resKOnonallelic <- as.data.table(results(ddsKO, contrast=c("allele","m","p"), altHypothesis="lessAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resKOnonallelic[padj<=pThresholds[1],]

# for genes imprinted in CTRL, test the difference:
dAx <- data.frame(genotype=gsub("^(.+)\\d_.+$","\\1",colnames(counts)), allele=gsub("^.+\\d_(.).+$","\\1",colnames(counts)), rep=gsub("^(.+\\d)_.+$","\\1",colnames(counts)), rep.n=as.factor(gsub("^.+(\\d)_.+$","\\1",colnames(counts))), stringsAsFactors=T)
rownames(dAx) <- colnames(counts)
jj <- rownames(counts)%in%resCTRLallelic[padj<=pThresholds[1], gene]
ddsX <- DESeqDataSetFromMatrix(countData = counts[jj,],
							  colData = dAx,
							  design = ~rep+allele) # this is what we would like but cannot do because the individual effect is nested in the groups: rep+genotype+allele+genotype:allele
design(ddsX) <- formula(~ rep.n + allele + genotype + genotype:rep.n + genotype:allele)
sizeFactors(ddsX) <- sizeFactors(ddsSample)[ceiling((1:ncol(counts))/2)]
ddsX <- DESeq(ddsX, betaPrior=F)	
resX <- as.data.table(results(ddsX, contrast=list("allelep.genotypematko", c("allele_p_vs_m","genotype_matko_vs_ctr")), altHypothesis="greaterAbs", lfcThreshold=0), keep.rownames="gene")
resX[padj<=pThresholds[1],]

# put all results together:
resChen <- merge(merge(resCTRLallelic, resKOallelic, by="gene", all=T, suffixes=c("_chen_ctrl", "_chen_mko")), resX, by="gene", all=T, suffixes=c("","_chen_x"))

X <- t(t(counts)/colSums(counts,na.rm=T))*1000000		
gs2 <- gs[gene%in%rownames(X),]	
pData1 <- data.table(
	gs2[,.(gene, list_type, log2fc_ours=log2fc)],
	ctr = log2(rowMeans(X[gs2$gene, grepl("ctr.+paternal",colnames(X))])+1) - log2(rowMeans(X[gs2$gene, grepl("ctr.+maternal",colnames(X))])+1),
	ko = log2(rowMeans(X[gs2$gene, grepl("matko.+paternal",colnames(X))])+1) - log2(rowMeans(X[gs2$gene, grepl("matko.+maternal",colnames(X))])+1)
)
tmp <- melt(as.data.table(counts, keep.rownames="gene"))
tmp[, sample_name:=gsub("_.+$","",variable)]
tmp[, allele:=ifelse(grepl("mat",gsub("^.+_","",variable), ignore.case=T),"m","p")]
tmp[, condition:=gsub("\\d+","",sample_name)]	
tmp2 <- dcast(tmp, condition+sample_name+gene ~ allele, value.var="value")
tmp2[, total:=m+p]
tmp2[, mat_ratio:=m/total]	

# add mKO data to main gene annotation table:
dtAll[rownames(counts), 	chen_ctr_p := rowSums(counts[, grepl("ctr.+paternal",colnames(counts))])]
dtAll[rownames(counts), 	chen_ctr_m := rowSums(counts[, grepl("ctr.+maternal",colnames(counts))])]
dtAll[rownames(counts), 	chen_ko_p := rowSums(counts[, grepl("matko.+paternal",colnames(counts))])]
dtAll[rownames(counts), 	chen_ko_m := rowSums(counts[, grepl("matko.+maternal",colnames(counts))])]
dtAll[rownames(X), 	chen_ctr := log2(rowMeans(X[, grepl("ctr.+paternal",colnames(X))])+1) - log2(rowMeans(X[, grepl("ctr.+maternal",colnames(X))])+1)]
dtAll[rownames(X), 	chen_ko := log2(rowMeans(X[, grepl("matko.+paternal",colnames(X))])+1) - log2(rowMeans(X[, grepl("matko.+maternal",colnames(X))])+1)]
dtAll[, chen_diff := chen_ko - chen_ctr]
colnames(counts) <- paste_("chen", colnames(counts))

dtAll <- merge(dtAll, as.data.table(counts, keep.rownames=T), by.x="gene", by.y="rn", all.x=T)





#--  Load second maternal KO dataset: --#


f <- resultsDir("Inoue_Supplemental_Table_S2.xlsx")
if(!file.exists(f)) {
	download.file("http://genesdev.cshlp.org/content/suppl/2018/11/21/gad.318675.118.DC1/Supplemental_Table_S2.xlsx", destfile=f)	
}
dExt <- as.data.table(readxl::read_xlsx(f, sheet=2, skip=1))[,-c(2:6)]
counts <- data.matrix(lib$dtToDf(dExt))
counts <- rbind(counts, "Snurf/Snrpn"=counts["Snurf",])
colnames(counts) <- tolower(colnames(counts))

dAx <- data.frame(genotype=gsub("^(.+)\\d_.+$","\\1",colnames(counts)), allele=gsub("^.+\\d_(.).+$","\\1",colnames(counts)), rep=gsub("^(.+\\d)_.+$","\\1",colnames(counts)), stringsAsFactors=F)
rownames(dAx) <- colnames(counts)

# estimate size factors on merged MAT and PAT counts:
sampleCounts <- counts[,c(1,3,5,7)] + counts[,c(2,4,6,8)]
sampleDA <- data.frame(genotype=gsub("^(.+)\\d_.+$","\\1",colnames(sampleCounts)), stringsAsFactors=F)
ddsSample <- DESeqDataSetFromMatrix(countData = sampleCounts,
							  colData = sampleDA,
							  design = ~genotype)
ddsSample <- DESeq(ddsSample, betaPrior=F)


j <- rownames(counts)%in%gs2$gene

# for CTRL:
i <- dAx$genotype=="ctr"
ddsCTRL <- DESeqDataSetFromMatrix(countData = counts[j,i],
							  colData = dAx[i,],
							  design = ~rep+allele)
sizeFactors(ddsCTRL) <- sizeFactors(ddsSample)[ceiling(which(i)/2)]
ddsCTRL <- DESeq(ddsCTRL, betaPrior=F)
resCTRLallelic <- as.data.table(results(ddsCTRL, contrast=c("allele","m","p"), altHypothesis="greaterAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resCTRLallelic[padj<=pThresholds[1],]
resCTRLnonallelic <- as.data.table(results(ddsCTRL, contrast=c("allele","m","p"), altHypothesis="lessAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resCTRLnonallelic[padj<=pThresholds[1],]

# for mKO:
i <- dAx$genotype!="ctr"
ddsKO <- DESeqDataSetFromMatrix(countData = counts[j,i],
							  colData = dAx[i,],
							  design = ~rep+allele)
sizeFactors(ddsKO) <- sizeFactors(ddsSample)[ceiling(which(i)/2)]
ddsKO <- DESeq(ddsKO, betaPrior=F)
resKOallelic <- as.data.table(results(ddsKO, contrast=c("allele","m","p"), altHypothesis="greaterAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resKOallelic[padj<=pThresholds[1],]
resKOnonallelic <- as.data.table(results(ddsKO, contrast=c("allele","m","p"), altHypothesis="lessAbs", independentFiltering=F, lfcThreshold=0.5), keep.rownames="gene")
resKOnonallelic[padj<=pThresholds[1],]

jj <- rownames(counts)%in%resCTRLallelic[padj<=pThresholds[1], gene]


# for genes imprinted in CTRL, test the difference:
dAx <- data.frame(genotype=gsub("^(.+)\\d_.+$","\\1",colnames(counts)), allele=gsub("^.+\\d_(.).+$","\\1",colnames(counts)), rep=gsub("^(.+\\d)_.+$","\\1",colnames(counts)), rep.n=as.factor(gsub("^.+(\\d)_.+$","\\1",colnames(counts))), stringsAsFactors=T)
rownames(dAx) <- colnames(counts)
ddsX <- DESeqDataSetFromMatrix(countData = counts[jj,],
							  colData = dAx,
							  design = ~rep+allele) # this is what we would like but cannot do because the individual effect is nested in the groups: rep+genotype+allele+genotype:allele
design(ddsX) <- formula(~ rep.n + allele + genotype + genotype:rep.n + genotype:allele)
sizeFactors(ddsX) <- sizeFactors(ddsSample)[ceiling((1:ncol(counts))/2)]
ddsX <- DESeq(ddsX, betaPrior=F)
resX <- as.data.table(results(ddsX, contrast=list("allelep.genotypematko", c("allele_p_vs_m","genotype_matko_vs_ctr")), altHypothesis="greaterAbs", lfcThreshold=0), keep.rownames="gene")
resX[padj<=pThresholds[1],]

resInoue <- merge(merge(resCTRLallelic, resKOallelic, by="gene", all=T, suffixes=c("_inoue_ctrl", "_inoue_mko")), resX, by="gene", all=T, suffixes=c("","_inoue_x"))

# put all results together:
resChenInoue <- merge(resChen, resInoue, all=T, suffixes=c("_chen_x","_inoue_x"))
fwrite(resChenInoue, file=resultsDir("chen_inoue_mko_deseq2.csv"))

X <- t(t(counts)/colSums(counts,na.rm=T))*1000000	
gs2 <- gs[gene%in%rownames(X),]	
pData2 <- data.table(
	gs2[,.(gene, list_type, log2fc_ours=log2fc)],
	ctr = log2(rowMeans(X[gs2$gene, grepl("CTR.+Pat",colnames(X))])+1) - log2(rowMeans(X[gs2$gene, grepl("CTR.+Mat",colnames(X))])+1),
	ko = log2(rowMeans(X[gs2$gene, grepl("matKO.+Pat",colnames(X))])+1) - log2(rowMeans(X[gs2$gene, grepl("matKO.+Mat",colnames(X))])+1)
)
tmp <- melt(as.data.table(counts, keep.rownames="gene"))
tmp[, sample_name:=gsub("_.+$","",variable)]
tmp[, allele:=ifelse(grepl("mat",gsub("^.+_","",variable), ignore.case=T),"m","p")]
tmp[, condition:=gsub("\\d+","",sample_name)]	
tmp2 <- dcast(tmp, condition+sample_name+gene ~ allele, value.var="value")
tmp2[, total:=m+p]
tmp2[, mat_ratio:=m/total]	
	
# add mKO data to main gene annotation table:
dtAll[rownames(counts), 	inoue_ctr_p := rowSums(counts[, grepl("CTR.+Pat",colnames(counts))])]
dtAll[rownames(counts), 	inoue_ctr_m := rowSums(counts[, grepl("CTR.+Mat",colnames(counts))])]
dtAll[rownames(counts), 	inoue_ko_p := rowSums(counts[, grepl("matKO.+Pat",colnames(counts))])]
dtAll[rownames(counts), 	inoue_ko_m := rowSums(counts[, grepl("matKO.+Mat",colnames(counts))])]
dtAll[rownames(X), 	inoue_ctr := log2(rowMeans(X[, grepl("CTR.+Pat",colnames(X), ignore.case=T)])+1) - log2(rowMeans(X[, grepl("CTR.+Mat",colnames(X), ignore.case=T)])+1)]
dtAll[rownames(X), 	inoue_ko := log2(rowMeans(X[, grepl("matKO.+Pat",colnames(X), ignore.case=T)])+1) - log2(rowMeans(X[, grepl("matKO.+Mat",colnames(X), ignore.case=T)])+1)]
dtAll[, inoue_diff := inoue_ko - inoue_ctr]

colnames(counts) <- paste_("inoue", colnames(counts))
dtAll <- merge(dtAll, as.data.table(counts, keep.rownames=T), by.x="gene", by.y="rn", all.x=T)
dtAll <- merge(dtAll, resChenInoue, by="gene", suffixes=c("","_2"), all=T)



# plot a compound heatmap of all results:

loadLibrary("ComplexHeatmap")
loadLibrary("circlize")

pData <- dtAll[gene%in%resChenInoue$gene,]

m <- lib$dtToDf(pData[, .(gene, chen_ctr, chen_ko, inoue_ctr, inoue_ko)]) # mean levels
m2 <- lib$dtToDf(pData[, .(gene, chen_diff, inoue_diff)]) # differences

qx <- ceiling(quantile(abs(m), 0.9, na.rm=T)) # cap at quantile for visualization
nBr <- 22
colFun <- colorRamp2(seq(-qx, qx, length=9), brewer.pal(9,"RdBu"))(seq(-qx, qx, length=21))

annot <- lib$dtToDf(gs[gene%in%rownames(m),.(gene,list_type)][order(list_type, gene),])
m <- m[rownames(annot),]

pvals <- data.matrix(lib$dtToDf(pData[, .(gene, padj_chen_ctrl, padj_chen_mko, padj_inoue_ctrl, padj_inoue_mko)]))[rownames(m),]
psigs <- t(apply(pvals, 1, lib$pToSig, ns="", threshs=pThresholds)) #n.s.
psigs[is.na(psigs)] <- ""
psigs[is.nan(pvals)] <- ""
pvals <- data.matrix(lib$dtToDf(pData[, .(gene, padj_chen_ctrl, padj_chen_mko, padj_inoue_ctrl, padj_inoue_mko)]))[rownames(m),]

aCols <- colorPalettes$genesetLabel[unique(annot[,1])]
aCols[["is_nbsx_not_nbix"]] <- alpha(aCols[["is_nbsx_not_nbix"]],0.25)

mx <- cbind(m, m2[rownames(m),])[,c(1,2,5,3,4,6)]
px <- cbind(psigs, naTo(pData[rownames(m),lib$pToSig(padj_chen_x,"n.s.")]), naTo(pData[rownames(m),lib$pToSig(padj_inoue_x,"n.s.")]))[,c(1,2,5,3,4,6)]
mx[px[,3]=="",3] <- NA
mx[px[,6]=="",6] <- NA
px2 <- px[px[,1]!="" | px[,4]!="" ,]
px2[px2=="n.s."] <- ""
annot2 <- annot[rownames(px2),,drop=F]
mx2 <- mx[rownames(px2),]

pheatmap::pheatmap(mx2, annotation_row=annot2, gaps_row=
as.numeric(cumsum(table(annot2[,1]))), display_numbers=px2, number_color="black",gaps_col=c(2,3,5), file=resultsDir("chen_inoue.pdf"), na_col = "#444444", col=colFun, breaks=seq(-qx, qx, length.out=nBr), border="white", cellwidth=16, cellheight=10, scale="none", cluster_cols=F, cluster_rows=F, annotation_colors=list(list_type=aCols))


pData <- cbind(annot2, mx2, pvals[rownames(annot2),], lib$dtToDf(pData[,.(gene,padj_chen_x, padj_inoue_x)])[rownames(annot2),])
fwrite(pData[pData$list_type!="is_unconfirmed",], file=resultsDir("source_data_figure_5a.csv"))
fwrite(pData, file=resultsDir("source_data_extended_data_figure_5a.csv"))

