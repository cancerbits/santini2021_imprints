#!/usr/bin/env Rscript
#
# Load "gold standard" GL-DMRs (list from Feil lab) and generate plots.
#
# run("imprinting", "gold_standard")
dependOn("imprinting", c("init"))

# Load somatic and germline DMRs:
feilSels <- sapply(c("use_gl_dmr","use_somatic_dmr"), function(feilSel) {

	dtFeil <- fread(dataDir("known_dmrs.csv"))[get(feilSel)==1,][,.(chrom,start,end, allele_imprinted, paternal, cpg_island, icr_name, name=make.unique(locus_name))]
	dtFeil[,k:=name]
	setkey(dtFeil,name)
	grFeil <- lib$dtToGr(dtFeil)

	regType <- paste_("feil",feilSel)

	# Calculate allelic coverage in regions:
	simpleCache(paste_("reads_and_meth", regType), {
		repeats <- NULL

		regSetX <- GRangesList(grFeil)
		names(regSetX) <- regType
		
		counts <- countAggregatedReads(dA, regSetX, name=paste_("rc",regType), repeats=repeats, sizes=sapply(regSetX,length))
		counts

	}, recreate=F, assignToVar="regReadCounts", reload=TRUE, cacheDir=analysisCacheDir)

	# Plot embryos/ESCs split into two heatmaps, but arranged in the same order:
	lib$pdfPlot(paste_("heatmap",regType,"custom"),13, 8)
	n <- "embryo"
	dWide <- makeWide(regReadCounts[id%in%sampleSelections[[n]] & id%in%dA$sample_name,], dtFeil)
	rowAnnot <- lib$dtToDf(dtFeil[rownames(dWide),.(ID=k, allele_imprinted, imprint_type=ifelse(paternal=="*","paternal","maternal"))][!is.na(ID),])
	colAnnot <- lib$dtToDf(dA[,.(sample_name, sample_group)])
	hm <- pheatmap(dWide*100, cellheight=10, cellwidth=10, breaks=seq(0,100,length.out=length(methCols)+1), show_rownames=T, cluster_rows=myHclust(dWide), cluster_cols=myHclust(t(dWide)), col=methCols, cutree_rows=2, treeheight_row=10, treeheight_col=10, annotation_row=rowAnnot[rownames(dWide),,drop=F], border_color="white", display_numbers=F, annotation_col=colAnnot, annotation_colors=colorPalettes[unique(c(colnames(rowAnnot), colnames(colAnnot)))])
	n <- "ESC"
	dWide <- makeWide(regReadCounts[id%in%sampleSelections[[n]] & id%in%dA$sample_name,], dtFeil)
	rowAnnot <- lib$dtToDf(dtFeil[rownames(dWide),.(ID=k, allele_imprinted, imprint_type=ifelse(paternal=="*","paternal","maternal"))][!is.na(ID),])
	colAnnot <- lib$dtToDf(dA[,.(sample_name, sample_group)])
	pheatmap(dWide*100, cellheight=10, cellwidth=10, breaks=seq(0,100,length.out=length(methCols)+1), show_rownames=T, cluster_rows=hm$tree_row, cluster_cols=myHclust(t(dWide)), col=methCols, cutree_rows=2, treeheight_row=0, treeheight_col=10, annotation_row=rowAnnot[rownames(dWide),,drop=F], border_color="white", display_numbers=F, annotation_col=colAnnot, annotation_colors=colorPalettes[unique(c(colnames(rowAnnot), colnames(colAnnot)))])
	dev.off()
	
	dWide <- makeWide(regReadCounts[id%in%unlist(sampleSelections[c("ESC","embryo")]) & id%in%dA$sample_name,], dtFeil)[hm$tree_row$order,]*100
	o <- dA[colnames(dWide),order(tissue,sample_group)]
	colnames(dWide) <- dA[colnames(dWide),][, sample_label_short]
	if(feilSel=="use_gl_dmr") {
		fwrite(dWide[,o], file=resultsDir("source_data_figure_2a.csv"))
	} else if(feilSel=="use_somatic_dmr") {
		fwrite(dWide[,o], file=resultsDir("source_data_extended_data_figure_2f.csv"))
	}
	
	feilOrder <- hm$tree_row$labels[hm$tree_row$order]
	
	list(dt=dtFeil, gr=grFeil, o=feilOrder)
}, simplify=F)

dtFeil <- feilSels[["use_gl_dmr"]]$dt
grFeil <- feilSels[["use_gl_dmr"]]$gr
feilOrder <- feilSels[["use_gl_dmr"]]$o
