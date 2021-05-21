#!/usr/bin/env Rscript
#
# Perform locus overlap analysis to annotate detected DMRs.
#
# run("imprinting", "region_enrichments")
dependOn("imprinting", c("define_regions","external","motifs"))

loadLibrary("LOLA")

# Define parameters:
effCol <- "log2odds"
filtBy <- "padj"
qThresh <- 0.005
effThresh <- log2(8)
tfCollections <- c("codex","encodeTFBSmm10","takahashi_2019")#"ucsc_features",

# Load annotation databases:
regionDB <- mergeRegionDBs(
	loadRegionDB(lib$resourceDir("regions/LOLACore/", genomeBuild), collections=tfCollections),
	loadRegionDB(lib$resourceDir("regions/chipseq_selected/", genomeBuild), collections=tfCollections) # we added BED files from Takahashi for ChIP-seq of Zfp445, https://doi.org/10.1101/gad.320069.118, available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115387
)
regionDB <- curateRegionDb(regionDB)

# Define user sets to test for enrichment:
lolaUserSets <- c(
	sapply(split(imprintedRegionCandidates[ranges_type!="dmrseq",], imprintedRegionCandidates[ranges_type!="dmrseq",ranges_type]), lib$dtToGr, "chr"),
	sapply(split(imprintedRegionCandidates[ranges_type=="dmrseq",], wangClustIds[imprintedRegionCandidates[ranges_type=="dmrseq",regionId]]), lib$dtToGr, "chr"),
	sapply(split(randomTrials, randomTrials$i), lib$dtToGr, "chr") 
)

selUserSets <- redefineUserSets(GRangesList(lolaUserSets) , regUniv)
nGrps <- sapply(selUserSets,length)

# Run LOLA (cache results to save time):
simpleCache("lolawith_custom_dmrseq_wang", {
	lolaResAll <- runLOLA(selUserSets, regUniv, regionDB, redefineUserSets=F)
	lolaResAll
}, assignToVar="lolaResAll", reload=T, recreate=F, cacheDir=analysisCacheDir)

lolaResAll[grepl("Wang",userSet), userSet:=wangClustNameMap[userSet]]

# Refine LOLA results by calculating p-values, log2 odds, tweaking term labels, etc:
lolaResAll <- lolaResAll[collection%in%tfCollections & userSet!="dmrseq",]
lolaResAll[, perc:=support/(support+c)]
lolaResAll[, percBg:=b/(b+d)]
lolaResAll[support==0, pValueLog:=0]
lolaResAll[, pval:=10^(-pValueLog)] # older versions of LOLA used log not log10 -- need to make sure this fits the version used (i.e. use exp(-p) for older versions instead)
lolaResAll[, padj:=p.adjust(pval, method="fdr")]
lolaResAll <- augmentLolaRes(lolaResAll, qThresh=qThresh, orCol="oddsRatio", qvalCol="padj", pvalCol="pval")
lolaResAll[, log2odds:=log2(oddsRatio)]
lolaResAll[support==0, log2odds:=lolaResAll[is.finite(log2odds), min(log2odds)]]
lolaResAll[b==0, log2odds:=lolaResAll[is.finite(log2odds), max(log2odds)]]

# Define significant enrichments:
lolaResAll[,sig:=get(filtBy)<=qThresh & get(effCol)>=effThresh]
lolaResAll[sig==T,.N,by=userSet]
sigI <-  lolaResAll[,which(sig & userSet%in%c(grep("DMR",names(nGrps),value=T)))] #,"feil"
selFiles <- lolaResAll[sigI,unique(filename)]
lolaResAll <- lolaResAll[filename%in%selFiles,]

# Select only terms enriched in ESC-datasets:
selTerms <- lolaResAll[!grepl("ref|rnd",userSet) & cellType=="Embryonic stem cell" & filename%in%selFiles & sig==T,][,.(s=absmax(10*log10(get(filtBy)))),by=term][order(s),term]  

selRowsFocus <-  lolaResAll[, unique(userSet)]
selTermsFocus <- selTerms

# Fix another mislabelling:
lolaResAll[cellType=="Embyonic stem cell" | cellType=="Es-e14", cellType:="Embryonic stem cell"]

# Plot results:
pData <- lolaResAll[term%in%selTerms,]
pData <- pData[filename%in%pData[support>0,unique(filename)],]
pData[, isESC:=ifelse(grepl("embryo",cellType,ignore.case=T),"ESC","Non-ESC")]
pData[, userSetX:=regionLabels[gsub("_[MP]$","",gsub("rnd_\\d+","rnd",userSet))]]
pData[is.na(userSetX), userSetX:=userSet]
pData[, nDatasets:=lib$dt2namedVec(pData[,length(unique(filename)),by=.(ID=term)])[term]]
pData[, minSites:=lib$dt2namedVec(pData[,min(size),by=.(ID=term)])[term]]
pData[, maxSites:=lib$dt2namedVec(pData[,max(size),by=.(ID=term)])[term]]
pData[, mxOr:=lib$dt2namedVec(pData[,max(log2odds),by=.(ID=term)])[term]]
pData[ , datasetNum:=as.numeric(as.factor(filename)), by=term]

selTermsFocus <-  pData[,max(mxOr),by=term][order(-V1),][1:7,term]

pDataF <- pData[term%in%selTermsFocus & filename%in%selFiles & isESC=="ESC",]

pDataF <- pDataF[, .(perc=max(perc), log2odds=max(log2odds), mxOr=max(mxOr)), by=.(term, userSet, userSetX, nDatasets)]

pLola <- ggplot(pDataF[!grepl("rnd_",userSet),], aes(x=reorder(lib$capFirst(tolower(term)),-perc), y=perc*100))  + geom_violin(data=pDataF[grepl("rnd_",userSet),], color=colorPalettes$regionLabel["Random regions"]) + geom_point(aes(shape=userSetX, color=userSetX)) + defTheme(noLegendTitle=T, topLegend=T) + coord_flip(clip="off") + ylab("Regions overlapping ChIP-seq peak (%)") + xlab(NULL) + scale_color_manual(values=colorPalettes$regionLabelWang) + scale_size_continuous(range=c(0.1, 4)) + scale_shape_manual(values=sapply(pDataF[,unique(userSetX)], function(x) ifelse(grepl("DMR-C",x),1,16)))

pdf(resultsDir("lola_bars_focus_wang.pdf"), 6, 3, useDingbats=F)
print(pLola)
dev.off()

fwrite(pDataF, file=resultsDir("source_data_figure_2e.csv"))



gg(patchwork::wrap_plots(pLoc + ggtitle("Genomic location"), pMotifs + ggtitle("Motif enrichments (HOCOMOCO v11, all significant)"), pLola + ggtitle("Overlaps with ChIP-seq peaks (LOLA, top 7)")) + patchwork::plot_layout(guides="collect", ncol=1) +  patchwork::plot_annotation(tag_levels = "A") & theme(legend.position = 'bottom'), "motifs_lola_pos", 7, 9, type="pdf")

# Write table with all results:
fwrite(pData[!grepl("rnd_",userSet),][order(-log2odds), .(userSet, collection, term, cellType, filename, a=support, b, c, d, log2odds, pval, padj)], file=resultsDir("lola_res_focus_wang.csv"))

