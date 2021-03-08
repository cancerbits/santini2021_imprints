#!/usr/bin/env Rscript
#
# Define clusters of nearby imprinted genes.
#
# run("imprinting", "imprint_control_regs")
dependOn("imprinting", c("annotate","expression"))

# Get all imprinted genes and sort them by genomic position
imprintGenes <- dtAll[is_bix | is_bsx | published_imprint,]
imprintGenes <- imprintGenes[order(chrom,start),]

# Check neighbouring genes pairwise to see if they are within the given distance:
winSize <- 250000
imprintGenes[, clustNum:=1:.N]
for(i in 2:nrow(imprintGenes)) {
	if(imprintGenes[i, chrom]==imprintGenes[i-1, chrom] & imprintGenes[i, start]<=imprintGenes[i-1, end+winSize]) {
		imprintGenes[i, clustNum:=imprintGenes[i-1, clustNum]]
	}
}

# Aggregate genes into clusters:
geneClusts <- imprintGenes[clustNum %in% imprintGenes[,.N,by=clustNum][N>1,clustNum],]
geneClusts[, clustNum:=as.numeric(as.factor(clustNum))]

# Annotate clusters:
nonCodRnas <- setdiff(transcriptAnnotationExpr[grepl("R_", id),unique(gene)],transcriptAnnotationExpr[grepl("M_", id),unique(gene)])

icrCandidates <- geneClusts[,.(
	genes=paste(gene,collapse="; "), 
	is_bsx=paste(intersect(gene, geneLists[list_type=="is_bsx",gene]),collapse="; "), 
	is_bix=paste(intersect(gene, geneLists[list_type=="is_bix",gene]),collapse="; "), 
	is_confirmed=paste(intersect(gene, geneLists[list_type=="is_confirmed",gene]),collapse="; "), 
	is_unconfirmed=paste(intersect(gene, geneLists[list_type=="is_unconfirmed",gene]),collapse="; "), 
	noncodernas=paste(intersect(gene, nonCodRnas),collapse="; "), 
	allelicrna_logfc_gt0=paste(intersect(gene, dtExprComplete[log2fc>0,gene]),collapse="; "), 
	allelicrna_logfc_lt0=paste(intersect(gene, dtExprComplete[log2fc<0,gene]),collapse="; "), 
	coord=sprintf("%s:%s-%s", unique(chrom), lib$prettyIntegers(min(start)), lib$prettyIntegers(max(end))), nearest_imprint=paste(unique(dmr_region_id),collapse="; "), nearest_dist=min(dmr_dist),
	h3k27me3_imprint=paste(intersect(gene, dtAll[h3k27me3_imprint==T,gene]),collapse="; "), 
	chrom=unique(chrom),
	start=min(start),
	end=max(end)
), by=.(cluster_number=clustNum)]
for(lt in geneLists[,unique(list_type)]) {
	icrCandidates[, paste_("overlap",lt):= sapply(genes, function(x) {
		paste(intersect(strsplit(x,"\\s*;\\s*")[[1]], geneLists[list_type==lt,gene]),collapse="; ")
	})]
}

# Write results:
fwrite(icrCandidates, file=resultsDir("icr_candidates_bix_bsx.csv"))


loadLibrary("ggstance")

i <- 1
p <- lapply(1:nrow(icrCandidates), function(i) {
	x <- icrCandidates[i,]
	
	win <- 50000
	xs <- x$start-win
	xe <- x$end+win
	
	reg <- GRanges(x$chrom, IRanges(xs,xe))
	
	# M = -1, P = 1
	
	rs <- imprintedRegionCandidates[chr==x$chrom & end>=xs & start<=xe & ranges_type%in%c("feil","dmrseq"), ]
	rs[, n:=1:.N]
	rs[, imprint:=ifelse(grepl("_M",grp), "-1","1")]
	
	gs <- dtAll[chrom==x$chrom & end>=xs & start<=xe, ]
	gs[, n:=1:.N]
	gs[, imprint:=ifelse(is_bsx|is_bix|published_imprint, as.character(sign(log2fc)),"0")]
	
	tads <- rbind(
		lib$grToDt(reduce(lib$dtToGr(dtTopol[chrom==x$chrom & end>=xs & start<=xe & `64CSE_G2`>=threshTopolScore, ])))[,.(chrom, start, end, imprint=1)], # G2 = pat
		lib$grToDt(reduce(lib$dtToGr(dtTopol[chrom==x$chrom & end>=xs & start<=xe & `64CSE_G1`>=threshTopolScore, ])))[,.(chrom, start, end, imprint=-1)] # G1 = mat
		
	)[!is.na(chrom),]
	
	
	tads[, n:=1:.N]
	
	h3k <- rbind(
		data.table(lib$grToDt(h3k27me3mat[overlapsAny(h3k27me3mat,reg)]),imprint="-1"),
		data.table(lib$grToDt(h3k27me3pat[overlapsAny(h3k27me3pat,reg)]),imprint="1")
	)[!is.na(chrom),]
	h3k[, n:=1:.N]
	
	pData <- rbind(
		rs[,.(start, end, imprint, lbl=regionId, type="dmr")],
		gs[,.(start, end, imprint, lbl=gene, type="gene")],
		h3k[,.(start, end, imprint, type="h3k27me3")],
		tads[,.(start, end, imprint, type="tad")],
	fill=T)[!is.na(`start`),]
	
	ypos <- c("tad","dmr","h3k27me3","gene")
	
	pData[, n:=1:.N, by=type]
	pData[, y:= as.double(sapply(type, function(x) which(x==ypos))- (n-1)/.N), by=type ]
	pData[start<xs, start:=xs]
	pData[end>xe, end:=xe]
	
	p <- ggplot(pData, mapping=aes(x=start,y=y,xend=end,yend=y)) + annotate(geom="rect", xmin=x$start, ymin=-0.1, xmax=x$end, ymax=length(ypos)+0.1, color="transparent", fill="orange", alpha=0.2) + theme_bw()  + xlab(x$chrom) + ylab(NULL) + ggtitle(sprintf("ICR %d (+/- %.1f kb),\ntotal length = %.1f kb",x$cluster_number, win/1000, (xe-xs)/1000)) + geom_segment(aes(color=imprint), size=2) + geom_text_repel(aes(label=lbl), direction="x", data=pData[!is.na(lbl),]) + scale_x_continuous(limits=c(xs,xe), oob=scales::squish, expand=c(0,0)) + scale_y_continuous(expand=c(0,0), oob=scales::squish, breaks=(1:length(ypos))-0.9, limits=c(-0.1,length(ypos)+0.1), labels=ypos) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(color="black")) + scale_color_manual(values=c("-1"="red","0"="grey","1"="blue"), guide=F)

	gg(p, paste_("icr",x$cluster_number), 5, 5, type="pdf")
	
	p + theme(legend.title=element_blank())
})
gg(cowplot::plot_grid(plotlist=p, ncol=6), paste_("icr", "all"), 18, 15, type="pdf")


# export figure source data:
pData <- rblapply(1:nrow(icrCandidates), function(i) {
	x <- icrCandidates[i,]
	
	win <- 50000
	xs <- x$start-win
	xe <- x$end+win
	
	reg <- GRanges(x$chrom, IRanges(xs,xe))
	
	# M = -1, P = 1
	
	rs <- imprintedRegionCandidates[chr==x$chrom & end>=xs & start<=xe & ranges_type%in%c("feil","dmrseq"), ]
	rs[, n:=1:.N]
	rs[, imprint:=ifelse(grepl("_M",grp), "-1","1")]
	
	gs <- dtAll[chrom==x$chrom & end>=xs & start<=xe, ]
	gs[, n:=1:.N]
	gs[, imprint:=ifelse(is_bsx|is_bix|published_imprint, as.character(sign(log2fc)),"0")]
	
	tads <- rbind(
		lib$grToDt(reduce(lib$dtToGr(dtTopol[chrom==x$chrom & end>=xs & start<=xe & `64CSE_G2`>=threshTopolScore, ])))[,.(chrom, start, end, imprint=1)], # G2 = pat
		lib$grToDt(reduce(lib$dtToGr(dtTopol[chrom==x$chrom & end>=xs & start<=xe & `64CSE_G1`>=threshTopolScore, ])))[,.(chrom, start, end, imprint=-1)] # G1 = mat
		
	)[!is.na(chrom),]
	
	
	tads[, n:=1:.N]
	
	h3k <- rbind(
		data.table(lib$grToDt(h3k27me3mat[overlapsAny(h3k27me3mat,reg)]),imprint="-1"),
		data.table(lib$grToDt(h3k27me3pat[overlapsAny(h3k27me3pat,reg)]),imprint="1")
	)[!is.na(chrom),]
	h3k[, n:=1:.N]
	
	pData <- rbind(
		rs[,.(start, end, imprint, lbl=regionId, type="dmr")],
		gs[,.(start, end, imprint, lbl=gene, type="gene")],
		h3k[,.(start, end, imprint, type="h3k27me3")],
		tads[,.(start, end, imprint, type="tad")],
	fill=T)[!is.na(`start`),]
	
	ypos <- c("tad","dmr","h3k27me3","gene")
	
	pData[, n:=1:.N, by=type]
	pData[, y:= as.double(sapply(type, function(x) which(x==ypos))- (n-1)/.N), by=type ]
	pData[start<xs, start:=xs]
	pData[end>xe, end:=xe]
	
	pData
})

confId <- pData[lbl%in%geneLists[list_type=="is_confirmed",gene], unique(id)]
nbsxId <- pData[lbl%in%geneLists[list_type=="is_nbsx",gene], unique(id)]
unconfId <- pData[lbl%in%geneLists[list_type=="is_unconfirmed",gene], unique(id)]

# assign cluster labels as used in the paper:
icrLabels <- c(
 "cluster 1" # 1
, "cluster 13" # 2
, "cluster 14" # 3
, "cluster 15" # 4
, "cluster 23" # 5
, "cluster 2" # 6
, "cluster 16" # 7
, "cluster 24" # 8
, "cluster 25" # 9
, "cluster 28" # 10
, "cluster 17" # 11
, "cluster 29" # 12
, "cluster 3" # 13
, "cluster 4" # 14
, "cluster 5" # 15
, "cluster 30" # 16
, "cluster 6" # 17
, "cluster 7" # 18
, "cluster 8" # 19
, "cluster 9" # 20
, "cluster 26" # 21
, "cluster 18" # 22
, "cluster 32" # 23
, "cluster 19" # 24
, "cluster 20" # 25
, "cluster 10" # 26
, "cluster 31" # 27
, "cluster 11" # 28
, "cluster 21" # 29
, "cluster 22" # 30
, "cluster 12" # 31
, "cluster 27" # 32
)
pData[, cluster_label:=icrLabels[id]]


pData <- pData[,.(id, cluster_label, type, y, lbl, imprint, start, end)]
fwrite(pData[id%in%intersect(confId,nbsxId),], file=resultsDir("source_data_figure_6a.csv"))
fwrite(pData[id%in%setdiff(nbsxId,confId),], file=resultsDir("source_data_figure_6b.csv"))
fwrite(pData[id%in%setdiff(confId,nbsxId),], file=resultsDir("source_data_extended_data_figure_6a.csv"))
fwrite(pData[id%in%setdiff(unconfId, c(confId, nbsxId)),], file=resultsDir("source_data_extended_data_figure_6b.csv"))



# 
# Plot chromosome ideograms with DMR, DEG positions etc:
if(doAllPlots) {
	loadLibrary("chromPlot")
	data(mm10_gap)
	binSize <- 1000000
	
	lib$pdfPlot("chromPlot_degs_icr", 10, 8)
	chromPlot(
		bin=binSize, maxSegs=1,
		segment=icrCandidates[,.(Chrom=chrom, Start=start, End=end)], colSegments="#00008b", 
		stat=icrCandidates[,.(Chrom=chrom, Start=start, End=end, Colors="darkgreen", ID=paste0("#",gsub("cluster ","",icrLabels[cluster_number])))], statCol="Value",
		statName="Value", statTyp="n", 
		stack=T, segLwd=5, 
		bands=mm10_cytoBandIdeo,
		yAxis=T, chrSide=c(-1,-1,1,1, 1, 1, -1, -1), noHist=T
	)
	dev.off()
	
	lib$pdfPlot("chromPlot_dmrs", 10, 8)
	chromPlot(
		annot1=lib$dtToGr(imprintedRegionCandidates[ranges_type=="dmrseq",],"chr"), colAnnot1="#daa520",
		annot2=lib$dtToGr(imprintedRegionCandidates[ranges_type=="feil",],"chr"), colAnnot2="#00008b",
		bin=binSize, maxSegs=1,
		stack=T, segLwd=5,
		bands=mm10_cytoBandIdeo,
		yAxis=T, chrSide=c(-1,-1,1,1, 1, 1, -1, -1), noHist=T
	)
	dev.off()
	
	lib$pdfPlot("chromPlot_dmrs_bg", 10, 8)
	chromPlot(
		annot1=lib$dtToGr(regionBg,"chr"), colAnnot1="#EEEEEE",
		bin=binSize, maxSegs=1,
		stack=T, segLwd=5,
		bands=mm10_cytoBandIdeo,
		yAxis=T, chrSide=c(-1,-1,1,1, 1, 1, -1, -1), noHist=T
	)
	dev.off()
		
	lib$pdfPlot("chromPlot_degs", 10, 8)
	chromPlot(
		annot1=lib$dtToGr(geneAnnotationExpr[gene%in%geneLists[list_type=="is_bsx",gene],]), colAnnot1="#daa520",
		annot2=lib$dtToGr(geneAnnotationExpr[gene%in%geneLists[list_type=="is_confirmed",gene],]), colAnnot2="#00008b",
		bin=binSize, maxSegs=1,	
		stack=T, segLwd=5,
		bands=mm10_cytoBandIdeo,
		yAxis=T, chrSide=c(1,1,1,1, 1, 1, -1, -1), noHist=T
	)
	dev.off()
	
	lib$pdfPlot("chromPlot_degs_bg", 10, 8)
	chromPlot(
		annot1=lib$dtToGr(geneAnnotationExpr), colAnnot1="#EEEEEE",
		bin=binSize, maxSegs=1,	
		stack=T, segLwd=5,
		bands=mm10_cytoBandIdeo,
		yAxis=T, chrSide=c(1,1,1,1, 1, 1, -1, -1), noHist=T
	)
	dev.off()

	lib$pdfPlot("chromPlot_h3k27me3_peaks", 10, 8)
	chromPlot(
			bin=binSize, maxSegs=1,	
			segment=dtAll[h3k27me3_imprint==T,.(Chrom=chrom, Start=start-binSize/5, End=end+binSize/5)], colSegments="#cb181d", 
			stack=T, segLwd=5,
			bands=mm10_cytoBandIdeo,
			yAxis=T, chrSide=c(1,1,1,1, 1, 1, -1, -1), noHist=T
		)
	dev.off()
	
	pData <- rbind(
		data.table(type="icr", icrCandidates[,.(chrom, start, end, label=icrLabels[cluster_number])]),
		data.table(type="dmr", imprintedRegionCandidates[ranges_type=="dmrseq",.(chrom=chr,start,end,label=regionId)]),
		data.table(type="known_dmr", imprintedRegionCandidates[ranges_type=="feil",.(chrom=chr,start,end,label=regionId)]),
		data.table(type="bsx", geneAnnotationExpr[gene%in%geneLists[list_type=="is_bsx",gene],.(chrom,start,end,label=gene)]),
		data.table(type="published_deg", geneAnnotationExpr[gene%in%geneLists[list_type=="published_imprint",gene],.(chrom,start,end,label=gene)]),
		data.table(type="h3k27me3",dtAll[h3k27me3_imprint==T,.(chrom, start=start-binSize/5, end=end+binSize/5,label=FALSE)]),
		data.table(type="dmr_background", regionBg[,.(chrom=chr,start,end,label=NA)]),
		data.table(type="deg_background", geneAnnotationExpr[,.(chrom,start,end,label=NA)])
	)
	
	fwrite(pData, file=resultsDir("source_data_figure_6c.csv"))
	
	
}