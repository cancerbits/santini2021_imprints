# Perform tissue-specific gene enrichment using the package TissueEnrich
# Author: Laura Santini

#load packages
library(TissueEnrich)
library(dplyr)

#load table annotating all results
tab_genes <- read.csv("tab_genes.csv",
                      header = TRUE, 
                      sep = ",")

#select gene lists
nBiX <- dplyr::filter(tab_genes, is_nbix==TRUE)
nBsX <- dplyr::filter(tab_genes, is_nbsx==TRUE)
Published_confirmed <- dplyr::filter(tab_genes, is_confirmed==TRUE)
Published_unconfirmed <- dplyr::filter(tab_genes, is_unconfirmed==TRUE)
Equivalent <- dplyr::filter(tab_genes, is_equivalent_0.1==TRUE)


#select top 105 equivalent genes
Equivalent_ordered <- dplyr::arrange(Equivalent, padj_equality_log2fc_threshold.1)
Equivalent_top106 <- dplyr::slice(Equivalent_ordered, 1:105)

#Convert lists of genes to character variables
nBiX_ch <- as.character(nBiX$gene_symbol)
nBsX_ch <- as.character(nBsX$gene_symbol)
Published_confirmed_ch <- as.character(Published_confirmed$gene_symbol)
Published_unconfirmed_ch <- as.character(Published_unconfirmed$gene_symbol)
Equivalent_top106_ch <- as.character(Equivalent_top106$gene_symbol)
All_detected_ch <- as.character(tab_genes$gene_symbol)

#Create  GeneSet objects (must contain:
#1. list of input genes as character vector
#2. organism 
#3. type of gene ID)

#Create GeneSet object for background genes using all detected genes in the RNA-Seq analysis
bg<-GeneSet(geneIds=All_detected_ch,organism="Mus Musculus",geneIdType=SymbolIdentifier())

#Create GeneSet objects for all list of interest
gs<-GeneSet(geneIds=nBiX_ch,organism="Mus Musculus",geneIdType=SymbolIdentifier())
gs2<-GeneSet(geneIds=nBsX_ch,organism="Mus Musculus",geneIdType=SymbolIdentifier())
gs8<-GeneSet(geneIds=Published_confirmed_ch,organism="Mus Musculus",geneIdType=SymbolIdentifier())
gs5<-GeneSet(geneIds=Published_unconfirmed_ch,organism="Mus Musculus",geneIdType=SymbolIdentifier())
gs7<-GeneSet(geneIds=Equivalent_top106_ch,organism="Mus Musculus",geneIdType=SymbolIdentifier())

#Run Tissue Enrichment 
output <- teEnrichment(inputGenes = gs,
                        rnaSeqDataset = 3, #3 is for "Mouse ENCODE",
                        tissueSpecificGeneType = 1, # All tissue specific genes
                        backgroundGenes = bg)

output2 <- teEnrichment(inputGenes = gs2,
                        rnaSeqDataset = 3, #3 is for "Mouse ENCODE",
                        tissueSpecificGeneType = 1, # All tissue specific genes
                        backgroundGenes = bg)

output8 <- teEnrichment(inputGenes = gs8,
                        rnaSeqDataset = 3, #3 is for "Mouse ENCODE",
                        tissueSpecificGeneType = 1, # All tissue specific genes
                        backgroundGenes = bg)

output5 <- teEnrichment(inputGenes = gs5,
                        rnaSeqDataset = 3, #3 is for "Mouse ENCODE",
                        tissueSpecificGeneType = 1, # All tissue specific genes
                        backgroundGenes = bg)

output7 <- teEnrichment(inputGenes = gs7,
                        rnaSeqDataset = 3, #3 is for "Mouse ENCODE",
                        tissueSpecificGeneType = 1, # All tissue specific genes
                        backgroundGenes = bg)

#Visualize results
seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
head(enrichmentOutput)

#Repeat for all output objects


#Export csv with enrichment results for all gene lists

write.csv(enrichmentOutput, 
          "enrichmentOutput_nBiX.csv", 
          row.names = TRUE)

write.csv(enrichmentOutput2, 
          "enrichmentOutput_nBsX.csv", 
          row.names = TRUE)

write.csv(enrichmentOutput8, 
          "enrichmentOutput_Published_confirmed.csv", 
          row.names = TRUE)

write.csv(enrichmentOutput5, 
          "enrichmentOutput_Published_uc.csv", 
          row.names = TRUE)

write.csv(enrichmentOutput7, 
          "enrichmentOutput_Equivalent_top106.csv", 
          row.names = TRUE)


##Plot combined enrichments for different gene lists in one barchart
library(ggplot2)
library(wesanderson)

#Load csv with previous results combined (with Excel)
TEnrich_combined <- read.csv("TEnrich_combined.csv",
                             header = TRUE, 
                             sep = ",")
theme_set(theme_classic())

pdf("TEnrich_final.pdf", width=9, height=3)
ggplot(data=TEnrich_combined, aes(x=Tissue, y=fold.change, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 5))+
  scale_x_discrete(limits=c("E14.5-Brain", "Cortex", "Olfactory bulb", 
                            "Cerebellum", "E14.5-Placenta", "E14.5-Heart" 
  )) #change order in the x axis

dev.off()

#Export csv with combined data as Source data
write.csv(TEnrich_combined, 
          "source_data_figure_1g.csv", 
          row.names = TRUE)
