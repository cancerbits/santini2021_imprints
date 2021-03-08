#!/usr/bin/env Rscript
#
# Perform DNA sequence motif enrichment analysis using FIMO.
#
# run("imprinting", "motifs")
dependOn("imprinting", c("gold_standard","define_regions","external"))

# Define reference database of known motifs:
motifDbs <- list(
	"hocomoco_full_mouse_v11"= list(meme=resourceDir("motifs/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme"), anno=resourceDir("motifs/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv"))
)
dir.create(resultsDir("fimo_out/"), showWarnings=F)


dtFeil[, regionId:=k]
dtUniv[, regionId:=k]
inputs <- list(fasta="univ", fimo="univ", dt=dtUniv, gr=regUniv)

# Export DNA sequences underlying all peaks in the "universe" and
# use FIMO to search them for matches to known motifs:

msg("Export sequences...")
n <- "univ"
seqFile <- resultsDir("fimo_out/",inputs$fasta,".fa")
if(!file.exists(seqFile)) {
	loadLibrary(genomeRef)
	seqs <- getSeq(get(genomeRef), inputs$gr)
	seqFile <- writeFASTA(seqs, seqFile)
} else {								
	msg("\t* sequence file exists -- skipped: ", seqFile)
}
		
for(motifDb in names(motifDbs)) {
	motifFile <- motifDbs[[motifDb]]$meme

	if(file.exists(motifFile)) {
		msg("Scan for motif occurrences in ", motifDb, " (", motifFile, ")...")
		
		f <- resultsDir("fimo_out/",inputs$fimo,"_", motifDb, "_t0.001.txt.gz")
		if(!file.exists(f)) {			
			f <- runFIMO(seqFile, outFile=f, motifFile=motifFile, inputIsFile=TRUE)
		} else {								
			msg("\t* output file exists -- skipped: ",f)
		}
		f
	}
	else {
		msg("ERROR: Unknown motif database: ", motifDb, " (", motifFile, ")...")			
	}
}


### Analyze results of motif search: ###

motifDbSel <- names(motifDbs)[1]
annoMotif <- fread(motifDbs[[motifDbSel]]$anno, key="Model")
selMotifDbs <- c(motifDbSel)
motifPThresh <- 1e-4  # default on FIMO website


# Load results and perform enrichment analysis. Cache the results to save time:
simpleCache("motif_enrich", {

	# Read FIMO results:
	fimo <- rblapply(selMotifDbs, function(motifDb) {
		msgF("reading %s", motifDb)
		motifFile <- motifDbs[[motifDb]]$meme
		resultsFile <- resultsDir("fimo_out/",inputs$fimo,"_", motifDb, "_t0.001.txt.gz")
		fimo <- data.table(readFIMOResult(resultsFile, motifFile=motifFile, motifPThresh=motifPThresh))	
		fimo[, regionId:=inputs$dt[as.numeric(substring(fimo$seq,4)), regionId]]
		fimo
	}, "motifDb")
	fimo[,motifName:=NULL]
	
	# Define region sets to test for enrichment:
	regsX <- rbind(
		regionBg[,.(chr, start, end, width, regionId, imprintType="bg", ranges_type="bg", grp="bg")],
		promos[,.(chr, start, end, width, regionId, imprintType="bg", ranges_type="bg", grp="bg")],
		dtFeil[,.(chr=chrom, start, end, width=end-start, regionId, imprintType="bg", ranges_type="bg", grp="bg")],
		imprintedRegionCandidates,
		randomTrials[,.(chr, start, end, width=end-start, regionId, imprintType="rnd", ranges_type="rnd", grp=i)],
	fill=T)
	regsX[, grp:=gsub("_[MP]","",grp)]
	regsX[, simpleRegId:=gsub("\\.\\d+$","",regionId)] ## don't distinguish multiple promoters of the same gene
	
	regsX[,nOverlapsFeil:=countOverlaps(lib$dtToGr(regsX,"chr"), grFeil)]
	
	# Restructure FIMO results by region sets:
	o <- findOverlaps(inputs$gr,lib$dtToGr(regsX,"chr"))		
	fimoX <- data.table(univId=inputs$dt[queryHits(o),regionId], regsX[subjectHits(o), .(regionId, simpleRegId, ranges_type, grp)])	
	fimoX <- merge(fimo[,.(univId=regionId, motifDb, motifId)], fimoX, by="univId", allow.cartesian=T, all.x=T)	
	fimoX <- fimoX[, .N, by=.(id=simpleRegId, motifDb, motifId, ranges_type, grp)]
	
	nGrp <- lib$dt2namedVec(regsX[!is.na(grp),length(unique(simpleRegId)),by=.(ID=grp)])	
	n <- nGrp[["bg"]]
	knownHits <- regsX[nOverlapsFeil>0, simpleRegId]

	# Summarize motif hits per group:
	motifEnrich <- fimoX[!is.na(grp), .( hitsMod=length(unique(id)), newHitsMod=length(setdiff(unique(id),knownHits)) ), by=.(grp, motifId)]
	motifEnrich[, motifName:=annoMotif[motifId, `Transcription factor`]]
	motifEnrich[, hitsBg:=lib$dt2namedVec(motifEnrich[grp=="bg", .(ID=motifId, hitsMod)])[motifId]]
	motifEnrich[, totalMod:=as.numeric(nGrp[as.character(grp)])]
	motifEnrich[, percMod:=hitsMod/totalMod]
	motifEnrich[, percBg:=hitsBg/n]
	motifEnrich[, percOnlyBg:=(hitsBg-hitsMod)/(n-totalMod)]
	motifEnrich <- motifEnrich[grp!="bg",]
		
	# Perform Fisher's test to determine significance of enrichments:
	motifEnrich[, A:=hitsMod]
	motifEnrich[, B:=totalMod-hitsMod]
	motifEnrich[, C:=hitsBg-hitsMod]
	motifEnrich[, D:=n-totalMod-hitsBg+hitsMod]
	if(nrow(motifEnrich[A<0|B<0|C<0|D<0,])>0) stop("sanity check failed: negative table entries")
	fishRes <- apply(data.matrix(as.matrix(motifEnrich[,.(A,B,C,D)])), 1, function(x) {
		m <- matrix(x,nrow=2,byrow=T)
		fishRes <- fisher.test(m, alternative="greater")
		c(fishRes$estimate,fishRes$p.value)
	})

	motifEnrich[, `:=`(c("oddsRatio","pval"), list(fishRes[1,], fishRes[2,]))]
	motifEnrich[, log2odds:=log2(oddsRatio)]	
	motifEnrich[, padj:=p.adjust(pval,method="fdr")]
	
	motifEnrich
}, assignToVar="motifEnrich", recreate=F, cacheDir=analysisCacheDir)
	

# Select enriched motifs:
motifsSel <- motifEnrich[!grepl("ref|rnd",grp) & padj<=0.005 & log2odds>=log2(4),][order(grp,-log2odds),]
motifsSelIds <- motifsSel[, unique(motifId)]
# ... exclude low-confidence motifs and motifs that are also strongly enriched in promoters (independent of imprinting):
motifsSelFocus <- grep(".[AB]$",setdiff(motifsSelIds, dcast(motifEnrich, motifId~gsub("_\\d+","",grp), value.var="log2odds", fun.aggregate=max)[dmrseq<(ref_promoter+log2(1.5)), motifId]), value=T)

# Write results table:
fwrite(motifEnrich[motifId%in%motifsSelIds,.(grp, motifId, motifName, percMod, percBg, pval, padj)], file=resultsDir("motif_enrich.csv"))

# Plot results:
pData <- motifEnrich[motifId%in%motifsSelFocus,]
pData[, grpX:=regionLabels[gsub("_[MP]$","",gsub("rnd_\\d+","rnd",grp))]]
pData[, mxOr:=lib$dt2namedVec(pData[,max(oddsRatio),by=.(ID=motifId)])[motifId]]
p <- ggplot(pData[grpX!="Random regions",], aes(x=reorder(motifId,mxOr), y=log2odds, color=grpX)) + geom_hline(yintercept=0, linetype="dashed") + geom_violin(data=pData[grpX=="Random regions",]) + geom_point() + defTheme(noLegendTitle=T) + coord_flip(clip="off") 
p <- p + scale_shape_manual(values=c(ESC=19,"Non-ESC"=1)) + ylab("Odds ratio, log2") + xlab(NULL) + scale_color_manual(values=colorPalettes$regionLabel) + scale_size_continuous(range=c(0.1, 4))
pdf(resultsDir("motif_bars.pdf"), 8, 1.7, useDingbats=F)
print(p)
dev.off()

fwrite(motifEnrich[motifId%in%motifsSelFocus & grp%in%c("dmrseq","feil"),.(grp, motifId, hitsMod, newHitsMod)], file=resultsDir("motif_hit_numbers.csv"))



####

# Load results and perform enrichment analysis. Cache the results to save time:
simpleCache("motif_enrich_wang", {

	# Read FIMO results:
	fimo <- rblapply(selMotifDbs, function(motifDb) {
		msgF("reading %s", motifDb)
		motifFile <- motifDbs[[motifDb]]$meme
		resultsFile <- resultsDir("fimo_out/",inputs$fimo,"_", motifDb, "_t0.001.txt.gz")
		fimo <- data.table(readFIMOResult(resultsFile, motifFile=motifFile, motifPThresh=motifPThresh))	
		fimo[, regionId:=inputs$dt[as.numeric(substring(fimo$seq,4)), regionId]]
		fimo
	}, "motifDb")
	fimo[,motifName:=NULL]
	
	# Define region sets to test for enrichment:
	regsX <- rbind(
		regionBg[,.(chr, start, end, width, regionId, imprintType="bg", ranges_type="bg", grp="bg")],
		promos[,.(chr, start, end, width, regionId, imprintType="bg", ranges_type="bg", grp="bg")],
		dtFeil[,.(chr=chrom, start, end, width=end-start, regionId, imprintType="bg", ranges_type="bg", grp="bg")],
		imprintedRegionCandidates,
		imprintedRegionCandidates[ranges_type=="dmrseq",.(chr, start, end, width, regionId, imprintType="wang", ranges_type="wang", grp=wangClustIds[regionId])],
		randomTrials[,.(chr, start, end, width=end-start, regionId, imprintType="rnd", ranges_type="rnd", grp=i)],
	fill=T)
	regsX[, grp:=gsub("_[MP]","",grp)]
	regsX[, simpleRegId:=gsub("\\.\\d+$","",regionId)] ## don't distinguish multiple promoters of the same gene
	
	regsX[,nOverlapsFeil:=countOverlaps(lib$dtToGr(regsX,"chr"), grFeil)]
	
	# Restructure FIMO results by region sets:
	o <- findOverlaps(inputs$gr,lib$dtToGr(regsX,"chr"))		
	fimoX <- data.table(univId=inputs$dt[queryHits(o),regionId], regsX[subjectHits(o), .(regionId, simpleRegId, ranges_type, grp)])	
	fimoX <- merge(fimo[,.(univId=regionId, motifDb, motifId)], fimoX, by="univId", allow.cartesian=T, all.x=T)	
	fimoX <- fimoX[, .N, by=.(id=simpleRegId, motifDb, motifId, ranges_type, grp)]
	
	nGrp <- lib$dt2namedVec(regsX[!is.na(grp),length(unique(simpleRegId)),by=.(ID=grp)])	
	n <- nGrp[["bg"]]
	knownHits <- regsX[nOverlapsFeil>0, simpleRegId]

	# Summarize motif hits per group:
	motifEnrich <- fimoX[!is.na(grp), .( hitsMod=length(unique(id)), newHitsMod=length(setdiff(unique(id),knownHits)) ), by=.(grp, motifId)]
	motifEnrich[, motifName:=annoMotif[motifId, `Transcription factor`]]
	motifEnrich[, hitsBg:=lib$dt2namedVec(motifEnrich[grp=="bg", .(ID=motifId, hitsMod)])[motifId]]
	motifEnrich[, totalMod:=as.numeric(nGrp[as.character(grp)])]
	motifEnrich[, percMod:=hitsMod/totalMod]
	motifEnrich[, percBg:=hitsBg/n]
	motifEnrich[, percOnlyBg:=(hitsBg-hitsMod)/(n-totalMod)]
	motifEnrich <- motifEnrich[grp!="bg",]
		
	# Perform Fisher's test to determine significance of enrichments:
	motifEnrich[, A:=hitsMod]
	motifEnrich[, B:=totalMod-hitsMod]
	motifEnrich[, C:=hitsBg-hitsMod]
	motifEnrich[, D:=n-totalMod-hitsBg+hitsMod]
	if(nrow(motifEnrich[A<0|B<0|C<0|D<0,])>0) stop("sanity check failed: negative table entries")
	fishRes <- apply(data.matrix(as.matrix(motifEnrich[,.(A,B,C,D)])), 1, function(x) {
		m <- matrix(x,nrow=2,byrow=T)
		fishRes <- fisher.test(m, alternative="greater")
		c(fishRes$estimate,fishRes$p.value)
	})

	motifEnrich[, `:=`(c("oddsRatio","pval"), list(fishRes[1,], fishRes[2,]))]
	motifEnrich[, log2odds:=log2(oddsRatio)]	
	motifEnrich[, padj:=p.adjust(pval,method="fdr")]
	
	motifEnrich
}, assignToVar="motifEnrich", recreate=F, cacheDir=analysisCacheDir)
	

motifEnrich[grepl("Wang",grp), grp:=wangClustNameMap[grp]]


# Select enriched motifs:
motifsSel <- motifEnrich[!grepl("ref|rnd",grp) & padj<=0.005 & log2odds>=log2(8),][order(grp,-log2odds),]
motifsSelIds <- motifsSel[, unique(motifId)]
# ... exclude low-confidence motifs and motifs that are also strongly enriched in promoters (independent of imprinting):
motifsSelFocus <- grep(".[AB]$",motifsSelIds, value=T)
# Write results table:
fwrite(motifEnrich[motifId%in%motifsSelIds,.(grp, motifId, motifName, percMod, percBg, pval, padj)], file=resultsDir("motif_enrich_wang.csv"))

# Plot results:
pData <- motifEnrich[motifId%in%motifsSelFocus,]
pData[, grpX:=regionLabels[gsub("_[MP]$","",gsub("rnd_\\d+","rnd",grp))]]
pData[is.na(grpX), grpX:=grp]
pData[, mxOr:=lib$dt2namedVec(pData[,max(oddsRatio),by=.(ID=motifId)])[motifId]]
pMotifs <- ggplot(pData[grpX!="Random regions",], aes(x=reorder(motifName,percMod), y=percMod*100))  + geom_violin(data=pData[grpX=="Random regions",], color=colorPalettes$regionLabel["Random regions"]) + geom_point(aes(shape=grpX, color=grpX)) + defTheme(noLegendTitle=T, topLegend=T) + coord_flip(clip="off") + scale_shape_manual(values=sapply(pData[,unique(grpX)], function(x) ifelse(grepl("DMR-C",x),1,16))) + ylab("Regions with motif (%)") + xlab(NULL) + scale_color_manual(values=colorPalettes$regionLabelWang)
pdf(resultsDir("motif_bars_wang.pdf"), 6, 3, useDingbats=F)
print(pMotifs)
dev.off()

fwrite(pData, file=resultsDir("source_data_figure_2f.csv"))