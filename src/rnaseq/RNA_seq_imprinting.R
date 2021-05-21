# Allele-specific RNA-seq analysis
# Author: Fabian Titz-Teixeira

require(ggplot2)
require(ggforce)
require(DESeq2)
require(org.Mm.eg.db)
require(gdata)
require(biomaRt)
require(openxlsx)

### get names of predicted genes from refseq ####
## file downloaded from (http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/ncbiRefSeq.txt.gz)
ncbiRef <- as.data.frame((read.table("ncbiRefSeq.txt", sep="\t", stringsAsFactors = FALSE)))
ncbiRef <- ncbiRef[,c(2,13)]
ncbiRef$predicted <- NA
for(d in 1:dim(ncbiRef)[1]){
  if(d%%10000==0){
    message(d)
  }
  if(substr(ncbiRef[d,1],1,1)=="X"){
    ncbiRef[d,3] <- 1
  }
  else{
    ncbiRef[d,3] <- 0
  }
}

predicted_genes <- unique(ncbiRef[which(ncbiRef$predicted==1),2])
predicted_genes <- predicted_genes[which(!predicted_genes%in%unique(ncbiRef[which(ncbiRef$predicted==0),2]))]

### create count table from allelome files ####

## read in counts from first config run
## samples 1 & 2 from forward and reverse cross
## middle part of the path depends on parameters in config file
sheet_1_2_1_2 <- read.table(file = "results/2020_03_26_B6_x_cast_1xE3_5_bla_1_B6_x_cast_1xE3_5_bla_2_cast_x_B6_1xE3_5_bla_1_cast_x_B6_1xE3_5_bla_2_anno_bigger100_260320_nonoverlap.bed_2_1/debug/read_count_per_SNP.txt",
                            header = TRUE)
sheet_1_2_1_2$position <- paste(sheet_1_2_1_2$chr, sheet_1_2_1_2$SNP_pos)
sheet1 <- sheet_1_2_1_2[,c(14,5,6,7,8,9,10,11,12,13)]
sheet1[][is.na(sheet1)] <- 0
colnames(sheet1) <- c("position", "gene","for_1_B6xCAST_B6_reads",
                      "for_1_B6xCAST_CAST_reads",
                      "for_2_B6xCAST_B6_reads",
                      "for_2_B6xCAST_CAST_reads",
                      "rev_1_CASTxB6_B6_reads",
                      "rev_1_CASTxB6_CAST_reads",
                      "rev_2_CASTxB6_B6_reads",
                      "rev_2_CASTxB6_CAST_reads")
counts1 <- aggregate(sheet1[,c(3:10)],list(sheet1$gene),sum)
colnames(counts1) <- c("gene","for_1_B6xCAST_B6_reads",
                       "for_1_B6xCAST_CAST_reads",
                       "for_2_B6xCAST_B6_reads",
                       "for_2_B6xCAST_CAST_reads",
                       "rev_1_CASTxB6_B6_reads",
                       "rev_1_CASTxB6_CAST_reads",
                       "rev_2_CASTxB6_B6_reads",
                       "rev_2_CASTxB6_CAST_reads")
xchr_1 <- list(unique(sheet_1_2_1_2[which(sheet_1_2_1_2$chr == "chrX"),5]))

## read in counts from second config run
## samples 3 & 1 from forward and 3 & 5 from reverse cross
sheet_3_1_3_5 <- read.table(file = "results/2020_03_26_B6_x_cast_1xE3_5_bla_3_B6_x_cast_1xE3_5_bla_1_cast_x_B6_1xE3_5_bla_3_cast_x_B6_1xE3_5_bla_5_anno_bigger100_260320_nonoverlap.bed_2_1/debug/read_count_per_SNP.txt",
                            header = TRUE)
sheet_3_1_3_5$position <- paste(sheet_3_1_3_5$chr, sheet_3_1_3_5$SNP_pos)
sheet2 <- sheet_3_1_3_5[,c(14,5,6,7,8,9,10,11,12,13)]
sheet2[][is.na(sheet2)] <- 0
colnames(sheet2) <- c("position", "gene","for_3_B6xCAST_B6_reads",
                      "for_3_B6xCAST_CAST_reads",
                      "for_1_B6xCAST_B6_reads",
                      "for_1_B6xCAST_CAST_reads",
                      "rev_3_CASTxB6_B6_reads",
                      "rev_3_CASTxB6_CAST_reads",
                      "rev_5_CASTxB6_B6_reads",
                      "rev_5_CASTxB6_CAST_reads")
counts2 <- aggregate(sheet2[,c(3,4,7,8,9,10)],list(sheet2$gene),sum)
colnames(counts2) <- c("gene","for_3_B6xCAST_B6_reads",
                       "for_3_B6xCAST_CAST_reads",
                       "rev_3_CASTxB6_B6_reads",
                       "rev_3_CASTxB6_CAST_reads",
                       "rev_5_CASTxB6_B6_reads",
                       "rev_5_CASTxB6_CAST_reads")
xchr_2 <- list(unique(sheet_3_1_3_5[which(sheet_3_1_3_5$chr == "chrX"),5]))

## read in counts from second config run
## samples 1 & 2 from forward and 7 & 5 from reverse cross
sheet_1_2_7_5 <- read.table(file = "results/2020_03_26_B6_x_cast_1xE3_5_bla_1_B6_x_cast_1xE3_5_bla_2_cast_x_B6_1xE3_5_bla_7_cast_x_B6_1xE3_5_bla_5_anno_bigger100_260320_nonoverlap.bed_2_1/debug/read_count_per_SNP.txt",
                            header = TRUE)
sheet_1_2_7_5$position <- paste(sheet_1_2_7_5$chr, sheet_1_2_7_5$SNP_pos)
sheet3 <- sheet_1_2_7_5[,c(14,5,6,7,8,9,10,11,12,13)]
sheet3[][is.na(sheet3)] <- 0
colnames(sheet3) <- c("position", "gene","for_1_B6xCAST_B6_reads",
                      "for_1_B6xCAST_CAST_reads",
                      "for_2_B6xCAST_B6_reads",
                      "for_2_B6xCAST_CAST_reads",
                      "rev_7_CASTxB6_B6_reads",
                      "rev_7_CASTxB6_CAST_reads",
                      "rev_5_CASTxB6_B6_reads",
                      "rev_5_CASTxB6_CAST_reads")
counts3 <- aggregate(sheet3[,c(7,8)],list(sheet3$gene),sum)
colnames(counts3) <- c("gene",
                       "rev_7_CASTxB6_B6_reads",
                       "rev_7_CASTxB6_CAST_reads")
xchr_3 <- list(unique(sheet_1_2_7_5[which(sheet_1_2_7_5$chr == "chrX"),5]))

## merge counts from all tables
merged_counts <- merge.data.frame(counts1, counts2, by="gene",all=TRUE)
merged_counts <- merge.data.frame(merged_counts, counts3, by="gene",all=TRUE)
merged_counts <- merged_counts[,c(1:5,10,11,6:9,12:17)]
merged_counts[][is.na(merged_counts)] <- 0
merged_counts <- data.frame(merged_counts[,-1], row.names = merged_counts[,1])

## subsetting the count table
counts_nox <- merged_counts[!row.names(merged_counts)%in%xchr_3[[1]],]
counts_nox <- counts_nox[!row.names(counts_nox)%in%xchr_2[[1]],]
counts_nox <- counts_nox[!row.names(counts_nox)%in%xchr_1[[1]],]


counts_mx <- counts_nox[apply(counts_nox[,], MARGIN = 1, function(x) any(x>=10)),]
counts_mx <- counts_mx[which(!row.names(counts_mx)%in%predicted_genes),]


### DESeq2 for imrpinting ####
design <- data.frame(row.names=c("for_1_B6xCAST_B6_reads",
                                 "for_1_B6xCAST_CAST_reads",
                                 "for_2_B6xCAST_B6_reads",
                                 "for_2_B6xCAST_CAST_reads",
                                 "for_3_B6xCAST_B6_reads",
                                 "for_3_B6xCAST_CAST_reads",
                                 "rev_1_CASTxB6_B6_reads",
                                 "rev_1_CASTxB6_CAST_reads",
                                 "rev_2_CASTxB6_B6_reads",
                                 "rev_2_CASTxB6_CAST_reads",
                                 "rev_3_CASTxB6_B6_reads",
                                 "rev_3_CASTxB6_CAST_reads",
                                 "rev_5_CASTxB6_B6_reads",
                                 "rev_5_CASTxB6_CAST_reads",
                                 "rev_7_CASTxB6_B6_reads",
                                 "rev_7_CASTxB6_CAST_reads"),
                     strain=c("B6","CAST","B6","CAST",
                              "B6","CAST","B6","CAST",
                              "B6","CAST","B6","CAST",
                              "B6","CAST","B6","CAST"),
                     imprint=c("maternal","paternal",
                               "maternal","paternal",
                               "maternal","paternal",
                               "paternal","maternal",
                               "paternal","maternal",
                               "paternal","maternal",
                               "paternal","maternal",
                               "paternal","maternal"),
                     sex=c("male","male","male",
                           "male","male","male",
                           "male","male","female",
                           "female","female","female",
                           "female","female","female",
                           "female"))


## maternal vs paternal
ddsFullCountTable_mx <- DESeqDataSetFromMatrix(countData = counts_mx,
                                               colData = design,
                                               design= ~ strain+imprint)
dds_mx <- DESeq(ddsFullCountTable_mx)
res_mx <- results(dds_mx, cooksCutoff = FALSE,name = "imprint_paternal_vs_maternal")
res_mx_ordered <- res_mx[order(res_mx$padj),]
res_mx_sign <- subset(res_mx_ordered, padj < 0.1)

## strain specific
ddsFullCountTable_strain <- DESeqDataSetFromMatrix(countData = counts_mx,
                                                   colData = design,
                                                   design= ~ imprint+strain)
dds_strain <- DESeq(ddsFullCountTable_strain)
res_strain <- results(dds_strain, cooksCutoff = FALSE)
res_strain_ordered <- res_strain[order(res_strain$padj),]

## create table for strain vs parental expression bias
strain_table <- as.data.frame(res_strain_ordered)
strain_table$name <- row.names(strain_table)
imprint_table <- as.data.frame(res_mx_ordered)
imprint_table$name <- row.names(imprint_table)
scatter_table <- merge.data.frame(strain_table,imprint_table, by="name", all = TRUE)
scatter_table <- scatter_table[,c(1,3,7,9,13)]
scatter_table[][is.na(scatter_table)] <- 1
scatter_table$col <- "grey"

for(d in 1:dim(scatter_table)[1]){
  if(scatter_table[d,3] < 0.1 && scatter_table[d,5] < 0.1){
    scatter_table[d,6] <- "forest green"
  }
  if(scatter_table[d,3] < 0.1 && scatter_table[d,5] >= 0.1){
    scatter_table[d,6] <- "red"
  }
  if(scatter_table[d,3] >= 0.1 && scatter_table[d,5] < 0.1){
    scatter_table[d,6] <- "blue"
  }
}

plot(scatter_table$log2FoldChange.x, 
     scatter_table$log2FoldChange.y,
     xlab="B6\t\t\t\t\t\t\t\t\tCAST", ylab="maternal\t\t\t\t\tpaternal",
     col=scatter_table$col)

### retrieve known imprints from wamidex and geneimprint and add info to results ####
## list of known imprinted genes from Schulz, R. et al., Epigenetics 2008 (doi: 10.4161/epi.3.2.5900) 
## https://atlas.genetics.kcl.ac.uk
## not provided in our upload
wamidex_imprints <- as.data.frame(read.csv("wamidex_ext.csv"))[,1]
## list of known imprinted genes from http://www.geneimprint.com/ 
## provided in our upload
geneimprint <- as.data.frame(read.csv("geneimprint_ext.csv"))[,1]

wiki_table <- as.data.frame(res_mx_ordered)
wiki_table$X <- row.names(wiki_table)
row.names(wiki_table) <- c()

wiki_table <- wiki_table[,c(7,1,2,6)]

wiki_table$cutoff20 <- 0
wiki_table$FPKM <- NA
wiki_table$inWAMIDEX <- 0
wiki_table$inGENEIMPRINT <- 0

## which genes do not have at least 20 reads covering SNPs in at least one sample
counts_mx_20 <- counts_nox[apply(counts_nox[,], MARGIN = 1, function(x) any(x>=20)),]

for(d in 1:dim(wiki_table)[1]){
  if(wiki_table[d,1]%in%row.names(counts_mx_20)){
    wiki_table[d,5] <- 1
  }
  if(wiki_table[d,1]%in%wamidex_imprints$Symbol | wiki_table[d,1]%in%wamidex_imprints$V3){
    wiki_table[d,7] <- 1
  }
  if(wiki_table[d,1]%in%geneimprint$Symbol | wiki_table[d,1]%in%geneimprint$V3){
    wiki_table[d,8] <- 1
  }
}

snp_table <- counts_mx[which(row.names(counts_mx)%in%wiki_table$X),]
snp_table$X <- row.names(snp_table)
wiki_table <- merge.data.frame(wiki_table,snp_table, by="X")


wiki_table[c(1:20),c(1,5,7,8)]
round(wiki_table[c(1:5),c(2:4)], 3)
colnames(wiki_table) <- c("gene", "baseMean", "log2fc", "padj", "cutoff20",
                          "FPKM", "in WAMIDEX", "in GENEIMPRINT",
                          "forward 1 B6 reads", "forward 1 CAST reads",
                          "forward 2 B6 reads", "forward 2 CAST reads",
                          "forward 3 B6 reads", "forward 3 CAST reads",
                          "reverse 1 B6 reads", "reverse 1 CAST reads",
                          "reverse 2 B6 reads", "reverse 2 CAST reads",
                          "reverse 3 B6 reads", "reverse 3 CAST reads",
                          "reverse 5 B6 reads", "reverse 5 CAST reads",
                          "reverse 7 B6 reads", "reverse 7 CAST reads")

### add ratio of contradicting snps and reads ####
## get counts per SNP for all SNPs found in at least one of the config runs
all_snp_counts <- merge.data.frame(sheet1[which(sheet1$gene%in%row.names(res_mx_ordered)),], 
                                   sheet2[which(sheet2$gene%in%row.names(res_mx_ordered)),], by="position",all=TRUE)
all_snp_counts <- merge.data.frame(all_snp_counts, sheet3[which(sheet3$gene%in%row.names(res_mx_ordered)),], 
                                   by="position",all=TRUE)

all_snp_counts[,2] <- as.character(all_snp_counts[,2])
for(d in 1:dim(all_snp_counts)[1]){
  if(d%%100000==0){
    message(d)
  }
  if(is.na(all_snp_counts[d,2])){
    if(is.na(all_snp_counts[d,11])){
      all_snp_counts[d,2] <- as.character(all_snp_counts[d,20])
    }
    else {
      all_snp_counts[d,2] <- as.character(all_snp_counts[d,11])
    }
  }
  if(all(is.na(all_snp_counts[d,c(2,11,20)]))){
    print(all_snp_counts[d,c(2,11,20)])
  }
}
all_snp_counts[,2] <- as.factor(all_snp_counts[,2])


all_snp_counts <- all_snp_counts[,c(1:6,12,13,7:10,16,17,18,19,25,26)]
colnames(all_snp_counts) <- c( "position","gene",                    
                               "for_1_B6xCAST_B6_reads","for_1_B6xCAST_CAST_reads",
                               "for_2_B6xCAST_B6_reads","for_2_B6xCAST_CAST_reads",
                               "for_3_B6xCAST_B6_reads","for_3_B6xCAST_CAST_reads",
                               "rev_1_CASTxB6_B6_reads","rev_1_CASTxB6_CAST_reads", 
                               "rev_2_CASTxB6_B6_reads","rev_2_CASTxB6_CAST_reads",  
                               "rev_3_CASTxB6_B6_reads","rev_3_CASTxB6_CAST_reads",  
                               "rev_5_CASTxB6_B6_reads","rev_5_CASTxB6_CAST_reads",
                               "rev_7_CASTxB6_B6_reads","rev_7_CASTxB6_CAST_reads")
all_snp_counts[][is.na(all_snp_counts)] <- 0
all_snp_counts[,c(3:dim(all_snp_counts)[2])] <- all_snp_counts[,c(3:dim(all_snp_counts)[2])]+1

## get counts on gene level
gene_counts_all <- wiki_table[,c(1,9:dim(wiki_table)[2])]
gene_counts_all[,c(2:dim(gene_counts_all)[2])] <- gene_counts_all[,c(2:dim(gene_counts_all)[2])]+1

head(gene_counts_all)

snp_plot_data_all <- all_snp_counts[,c(1,2)]
snp_plot_data_all$snpratio <- 0
snp_plot_data_all$generatio <- 0
snp_plot_data_all$meansnp <- 0
snp_plot_data_all$meangene <- 0
snp_plot_data_all$col <- "black"

## get rations per SNP and per gene (paternal vs maternal)
for(d in 1:dim(snp_plot_data_all)[1]){
  snp_plot_data_all[d,3] <- log2(mean(c(rowMeans(all_snp_counts[d,c(4,6,8)])/
                                          rowMeans(all_snp_counts[d,c(3,5,7)]),
                                        rowMeans(all_snp_counts[d,c(9,11,13,15,17)])/
                                          rowMeans(all_snp_counts[d,c(10,12,14,16,18)]))))
  snp_plot_data_all[d,4] <- log2(mean(c(rowMeans(gene_counts_all[which(gene_counts_all$gene == as.character(snp_plot_data_all[d,2])),c(3,5,7)])/
                                          rowMeans(gene_counts_all[which(gene_counts_all$gene == as.character(snp_plot_data_all[d,2])),c(2,4,6)]),
                                        rowMeans(gene_counts_all[which(gene_counts_all$gene == as.character(snp_plot_data_all[d,2])),c(8,10,12,14,16)])/
                                          rowMeans(gene_counts_all[which(gene_counts_all$gene == as.character(snp_plot_data_all[d,2])),c(9,11,13,15,17)]))))
  snp_plot_data_all[d,5] <- mean(rowMeans(all_snp_counts[d,c(3:18)]))
  snp_plot_data_all[d,6] <- mean(rowMeans(gene_counts_all[which(gene_counts_all$gene == as.character(snp_plot_data_all[d,2])),c(2:17)]))
  if(scatter_table[which(scatter_table$name == as.character(snp_plot_data_all[d,2])),6] == "forest green"){
    snp_plot_data_all[d,7] <- "forest green"
  }
  
}

## which SNPs contradict the ratio on the gene level
not_trusted_all <- snp_plot_data_all
not_trusted_all$contradict <- 0
for(d in 1:dim(not_trusted_all)[1]){
  if(not_trusted_all[d,4] > 0){
    if(not_trusted_all[d,3] < 0){
      not_trusted_all[d,8] <- 1
    }
  }
  if(not_trusted_all[d,4] < 0){
    if(not_trusted_all[d,3] > 0){
      not_trusted_all[d,8] <- 1
    }
  }
  if(not_trusted_all[d,4] == 0){
    if(not_trusted_all[d,3] != 0){
      not_trusted_all[d,8] <- 1
    }
  }
}

contradiction_ratios_all <- data.frame(unique(not_trusted_all$gene))
contradiction_ratios_all$ratio <- 0
contradiction_ratios_all$readratio <- 0
colnames(contradiction_ratios_all) <- c("gene", "ratio", "readratio")

## get ratios
for(d in 1:dim(contradiction_ratios_all)[1]){
  contr_sum <- 0
  sup_sum <- 0
  contradiction_ratios_all[d,2] <- sum(not_trusted_all[which(not_trusted_all$gene == contradiction_ratios_all$gene[d]),8])/dim(not_trusted_all[which(not_trusted_all$gene == contradiction_ratios_all$gene[d]),])[1]
  
  data_contr <- not_trusted_all[which(not_trusted_all$gene == contradiction_ratios_all$gene[d]),]
  for(r in 1:dim(data_contr)[1]){
    if(data_contr[r,8] == 1){
      contr_sum <- contr_sum + data_contr[r,5]-1
    }
    else{
      sup_sum <- sup_sum + data_contr[r,5]-1
    }
  }
  contradiction_ratios_all[d,3] <- contr_sum/sum(contr_sum,sup_sum)
}

wiki_final_imprint <- merge.data.frame(wiki_table, contradiction_ratios_all, by = "gene", all = FALSE)


### add statistics to test for equevalence ####
dds_mx_eq <- DESeq(ddsFullCountTable_mx, betaPrior = FALSE)
res_mx_eq_1 <- results(dds_mx_eq, cooksCutoff = FALSE, lfcThreshold = 1, altHypothesis = "lessAbs")
res_mx_ordered_eq_1 <- res_mx_eq_1[order(res_mx_eq_1$padj),]
res_mx_sign_eq_1 <- subset(res_mx_ordered_eq_1, padj < 0.1)
res_mx_sign_eq_1

res_mx_1 <- results(dds_mx_eq, cooksCutoff = FALSE, lfcThreshold = 1, altHypothesis = "greaterAbs")
res_mx_ordered_1 <- res_mx_1[order(res_mx_1$padj),]
res_mx_sign_1 <- subset(res_mx_ordered_1, padj < 0.1)
res_mx_sign_1

df_mx_ordered_1 <- as.data.frame(res_mx_ordered_1)
df_mx_ordered_1$gene <- row.names(df_mx_ordered_1)
row.names(df_mx_ordered_1) <- c()
df_mx_ordered_1 <- df_mx_ordered_1[,c(7,6)]
colnames(df_mx_ordered_1) <- c("gene", "padj differential lg2fcThreshold 1")

df_mx_ordered_eq_1 <- as.data.frame(res_mx_ordered_eq_1)
df_mx_ordered_eq_1$gene <- row.names(df_mx_ordered_eq_1)
row.names(df_mx_ordered_eq_1) <- c()
df_mx_ordered_eq_1 <- df_mx_ordered_eq_1[,c(7,6)]
colnames(df_mx_ordered_eq_1) <- c("gene", "padj equality lg2fcThreshold 1")

wiki_final_imprint <- merge.data.frame(wiki_final_imprint, df_mx_ordered_1, by = "gene")
wiki_final_imprint <- merge.data.frame(wiki_final_imprint, df_mx_ordered_eq_1, by = "gene")


wiki_final_imprint <- wiki_final_imprint[,c(1:5,7,8,25:28,9:24)]


colnames(wiki_final_imprint) <- c("gene", "baseMean", "log2fc", "padj", "cutoff20",
                                  "in WAMIDEX", "in GENEIMPRINT", 
                                  "ratio of contradicting SNPs",
                                  "ratio of contradicting reads",
                                  "padj differential lg2fcThreshold 1",
                                  "padj equality lg2fcThreshold 1",
                                  "forward 1 B6 maternal reads", "forward 1 CAST paternal reads",
                                  "forward 2 B6 maternal reads", "forward 2 CAST paternal reads",
                                  "forward 3 B6 maternal reads", "forward 3 CAST paternal reads",
                                  "reverse 1 B6 paternal reads", "reverse 1 CAST maternal reads",
                                  "reverse 2 B6 paternal reads", "reverse 2 CAST maternal reads",
                                  "reverse 3 B6 paternal reads", "reverse 3 CAST maternal reads",
                                  "reverse 5 B6 paternal reads", "reverse 5 CAST maternal reads",
                                  "reverse 7 B6 paternal reads", "reverse 7 CAST maternal reads")

wiki_sign <- wiki_final_imprint[which(wiki_final_imprint$padj<0.1),]

## which genes do not have at least 20 reads in all samples from one cross
reject <- c()
for(d in 1:dim(wiki_sign)[1]){
  if(sum(wiki_sign[d,c(12:17)])<20){
    reject <- c(reject, as.character(wiki_sign[d,1]))
  }
  if(sum(wiki_sign[d,c(18:27)])<20){
    reject <- c(reject, as.character(wiki_sign$gene[d]))
  }
}

## which genes have conflicting direction of imprinting between forward and reverse cross
contra <- c()
for(d in 1:dim(wiki_final_imprint)[1]){
  forw <- log2((sum(wiki_final_imprint[d,c(12,14,16)])+1)/(sum(wiki_final_imprint[d,c(13,15,17)])+1))
  rev <- log2((sum(wiki_final_imprint[d,c(19,21,23,25,27)])+1)/(sum(wiki_final_imprint[d,c(18,20,22,24,26)])+1))
  
  if(forw == 0|rev == 0){
    contra <- c(contra,as.character(wiki_final_imprint[d,1]))
  }
  else{
    if(forw < 0 & rev > 0){
      contra <- c(contra,as.character(wiki_final_imprint[d,1]))
    }
    if(forw > 0 & rev < 0){
      contra <- c(contra,as.character(wiki_final_imprint[d,1]))
    }
  }
}

reject <- unique(c(reject,contra))

wiki_final_imprint$reject <- 0
wiki_final_imprint[which(wiki_final_imprint$gene%in%reject),28] <- 1
wiki_final_imprint <- wiki_final_imprint[,c(1:9, 28, 10:27)]

### 147 identified parent-of-origin-specific genes
dim(wiki_final_imprint[which(wiki_final_imprint$padj<0.1&
                               wiki_final_imprint$`ratio of contradicting reads`<0.2&
                               wiki_final_imprint$reject==0),])

### add imprints from Otago to outputs ####
## list of known imprinted genes from www.otago.ac.nz/IGC
## not provided in our upload
Otago_imprints <- as.data.frame(read.csv("Otago.csv",header = FALSE))

egSYMBOL <- toTable(org.Mm.egALIAS2EG)


for(d in 1:dim(Otago_imprints)[1]){
  if (dim(egSYMBOL[which(egSYMBOL$alias_symbol == as.character(Otago_imprints[d,1])),])[1]==1) {
    Otago_imprints[d,2] <- get(as.character(Otago_imprints[d,1]), org.Mm.egALIAS2EG)
    Otago_imprints[d,3] <- get(as.character(Otago_imprints[d,2]), org.Mm.egSYMBOL)
  }
}

Otago_imprints[which(Otago_imprints$V3 == ""),2] <- NA

wiki_final_imprint$in.Otago <- 0

for(d in 1:dim(wiki_final_imprint)[1]){
  if(wiki_final_imprint[d,1]%in%Otago_imprints$V1 | wiki_final_imprint[d,1]%in%Otago_imprints$V3){
    wiki_final_imprint[d,29] <- 1
  }
}



wiki_final_imprint <- wiki_final_imprint[,c(1:7,29,8:28)]

### adding manually curated information to imprinting table ####
table_laura <- read.xls("Table_S1 imprinted genes list_LS_090320.xlsx")
table_laura <- table_laura[,1:15]


imprint_finished_table <- wiki_final_imprint

imprint_finished_table <- imprint_finished_table[,c(1:5,9:29)]

imprint_finished_table <- merge(imprint_finished_table, table_laura,
                                by.x="gene",by.y="Gene.symbol",all.x=T,all.y=F)

imprint_finished_table$cutoff20 <- as.numeric(apply(imprint_finished_table[,11:26],1,max)>=20)

#write.csv(imprint_finished_table,"wiki_final_imprints_full_results.csv",row.names = F)


### finalize tables for known imprints (different names in some lists) ####
## list of known imprinted genes from Schulz, R. et al., Epigenetics 2008 (doi: 10.4161/epi.3.2.5900) 
## https://atlas.genetics.kcl.ac.uk
## not provided in our upload
wamidex <- read.csv("wamidex_ext.csv")
wamidex[which(wamidex$Symbol == ""),2] <- NA
## list of known imprinted genes from http://www.geneimprint.com/ 
## provided in our upload
imprintlist <- read.csv("geneimprint_ext.csv")
imprintlist[which(imprintlist$Symbol == ""),2] <- NA

merged_imprint <- merge.data.frame(wamidex, imprintlist, 
                                   by.x = "Name", by.y = "name", all = TRUE)
## list of known imprinted genes from Blake, A. et al., Nucleic Acids 2010 (doi:10.1093/nar/gkp867)
## https://www.mousebook.org/
## not provided in our upload
mousebook <- read.xls("mousebook imprinted genes.xlsx")
mousebook[which(mousebook$Gene == ""),2] <- NA


merged_imprint$Symbol <- NA
for(d in 1:dim(merged_imprint)[1]){
  if(is.na(merged_imprint[d,2])){
    if(is.na(merged_imprint[d,3])){
      merged_imprint[d,4] <- NA
    }
    else{
      merged_imprint[d,4] <- as.character(merged_imprint[d,3])
    }
  }
  else{
    if(is.na(merged_imprint[d,3])){
      merged_imprint[d,4] <- as.character(merged_imprint[d,2])
    }
    else{
      merged_imprint[d,4] <- as.character(merged_imprint[d,2])
      #print(as.character(merged_imprint[3,2]) == as.character(merged_imprint[3,3]))
    }
  }
}

merged_imprint <- merged_imprint[,c(1,4)]

merged_imprint$Name <- as.character(merged_imprint$Name)
merged_imprint$mousebook <- ""

for(d in 1:dim(merged_imprint)[1]){
  if(merged_imprint[d,1]%in%mousebook$Gene|merged_imprint[d,2]%in%mousebook$Gene){
    merged_imprint[d,3] <- as.character(mousebook[which(mousebook$Gene==merged_imprint[d,1]|
                                                          mousebook$Gene==merged_imprint[d,2]),1])
  }
  else{
    merged_imprint[d,3] <- NA
  }
}

imprinted_list <- rbind(merged_imprint,data.frame(Name=NA,Symbol=NA,mousebook=mousebook$Gene[which(!mousebook$Gene%in%merged_imprint$mousebook)]))

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("entrezgene_id", "gene_biotype","mgi_symbol"),mart=ensembl)

mousebook_symbols <- egSYMBOL[which(egSYMBOL$alias_symbol%in%na.omit(imprinted_list$mousebook)),]
mousebook_lookup <- biomart[which(biomart$entrezgene_id%in%mousebook_symbols$gene_id),]
mousebook_mgi <- merge.data.frame(mousebook_lookup, mousebook_symbols, 
                                  by.x = "entrezgene_id", by.y = "gene_id", all.x = TRUE)

imprinted_list_fin <- imprinted_list
imprinted_list_fin[which(is.na(imprinted_list_fin$Name)),1] <- imprinted_list_fin[which(is.na(imprinted_list_fin$Name)),3]

for(d in which(is.na(imprinted_list$Name))){
  if(imprinted_list_fin[d,1]%in%mousebook_mgi[,4]){
    imprinted_list_fin[d,2] <- mousebook_mgi[which(mousebook_mgi$alias_symbol==imprinted_list_fin[d,1]),3][1]
  }
}

Otago_imprints$V1 <- as.character(Otago_imprints$V1)
Otago_imprints$V3 <- as.character(Otago_imprints$V3)


imprinted_list_fin$Otago <- NA

for(d in 1:dim(imprinted_list_fin)[1]){
  if(!is.na(imprinted_list_fin[d,2])){
    if(imprinted_list_fin[d,2]%in%Otago_imprints$V3){
      imprinted_list_fin[d,4] <- Otago_imprints[which(Otago_imprints$V3==imprinted_list_fin[d,2]),1][1]
    }
    if(imprinted_list_fin[d,1]%in%Otago_imprints$V1){
      imprinted_list_fin[d,4] <- Otago_imprints[which(Otago_imprints$V1==imprinted_list_fin[d,1]),1][1]
    }
  }
  else{
    if(imprinted_list_fin[d,1]%in%Otago_imprints$V1){
      imprinted_list_fin[d,4] <- Otago_imprints[which(Otago_imprints$V1==imprinted_list_fin[d,1]),1][1]
    }
  }
}

imprinted_list_fin <- rbind(imprinted_list_fin,data.frame(Name=Otago_imprints[which(!Otago_imprints$V1%in%imprinted_list_fin$Otago),1],
                                                                        Symbol=Otago_imprints[which(!Otago_imprints$V1%in%imprinted_list_fin$Otago),2],
                                                                        mousebook=NA,
                                                                        Otago=Otago_imprints[which(!Otago_imprints$V1%in%imprinted_list_fin$Otago),1]))

imprinted_list_fin <- imprinted_list_fin[order(imprinted_list_fin$Name),]


known_impr <- imprinted_list_fin[,1:2]

### merge all known imprinted genes from different lists into final table ####
head(wamidex) ##Wamidex
head(imprintlist) ##Geneimprint
head(Otago_imprints) ##Otago
head(impr_lit) ##Mousebook

final_literature_imprints <- known_impr

final_literature_imprints$Wamidex <- NA
final_literature_imprints$Geneimprint <- NA
final_literature_imprints$Otago <- NA
final_literature_imprints$Mousebook <- NA


final_literature_imprints[which(final_literature_imprints$Name%in%wamidex$Name),3] <- 1
final_literature_imprints[which(final_literature_imprints$Name%in%na.omit(wamidex$Symbol)),3] <- 1
final_literature_imprints[which(final_literature_imprints$Symbol%in%wamidex$Name),3] <- 1
final_literature_imprints[which(final_literature_imprints$Symbol%in%na.omit(wamidex$Symbol)),3] <- 1

final_literature_imprints[which(final_literature_imprints$Name%in%imprintlist$name),4] <- 1
final_literature_imprints[which(final_literature_imprints$Name%in%na.omit(imprintlist$Symbol)),4] <- 1
final_literature_imprints[which(final_literature_imprints$Symbol%in%imprintlist$name),4] <- 1
final_literature_imprints[which(final_literature_imprints$Symbol%in%na.omit(imprintlist$Symbol)),4] <- 1

final_literature_imprints[which(final_literature_imprints$Name%in%Otago_imprints$V1),5] <- 1
final_literature_imprints[which(final_literature_imprints$Name%in%na.omit(Otago_imprints$V3)),5] <- 1
final_literature_imprints[which(final_literature_imprints$Symbol%in%Otago_imprints$V1),5] <- 1
final_literature_imprints[which(final_literature_imprints$Symbol%in%na.omit(Otago_imprints$V3)),5] <- 1

final_literature_imprints[which(final_literature_imprints$Name%in%impr_lit$Gene),6] <- 1
final_literature_imprints[which(final_literature_imprints$Symbol%in%impr_lit$Gene),6] <- 1


final_literature_imprints[which(is.na(final_literature_imprints$Wamidex)),3] <- 0
final_literature_imprints[which(is.na(final_literature_imprints$Geneimprint)),4] <- 0
final_literature_imprints[which(is.na(final_literature_imprints$Otago)),5] <- 0
final_literature_imprints[which(is.na(final_literature_imprints$Mousebook)),6] <- 0


final_literature_imprints$Allele_from_mousebook <- NA
final_literature_imprints$Chr_from_mousebook <- NA
for(d in 1:dim(final_literature_imprints)[1]){
  if(final_literature_imprints[d,6]==1){
    final_literature_imprints[d,7] <- impr_lit[which(impr_lit$Gene%in%na.omit(c(as.character(final_literature_imprints[d,1]),
                                                                                as.character(final_literature_imprints[d,2])))),5]
    final_literature_imprints[d,8] <- impr_lit[which(impr_lit$Gene%in%na.omit(c(as.character(final_literature_imprints[d,1]),
                                                                                as.character(final_literature_imprints[d,2])))),4]
  }
}

final_literature_imprints[which(final_literature_imprints$Allele_from_mousebook==1),7] <- "M"
final_literature_imprints[which(final_literature_imprints$Allele_from_mousebook==2),7] <- "P"

imprints_mgi_only <- data.frame(mgi=na.omit(unique(final_literature_imprints$Symbol)),
                                Wamidex=0,Geneimprint=0,Otago=0,Mousebook=0)

for(d in 1:dim(imprints_mgi_only)[1]){
  if(colSums(final_literature_imprints[which(final_literature_imprints$Symbol==imprints_mgi_only$mgi[d]),3:6])[1]>0){
    imprints_mgi_only[d,2] <- 1
  }
  if(colSums(final_literature_imprints[which(final_literature_imprints$Symbol==imprints_mgi_only$mgi[d]),3:6])[2]>0){
    imprints_mgi_only[d,3] <- 1
  }
  if(colSums(final_literature_imprints[which(final_literature_imprints$Symbol==imprints_mgi_only$mgi[d]),3:6])[3]>0){
    imprints_mgi_only[d,4] <- 1
  }
  if(colSums(final_literature_imprints[which(final_literature_imprints$Symbol==imprints_mgi_only$mgi[d]),3:6])[4]>0){
    imprints_mgi_only[d,5] <- 1
  }
}

#write.csv(imprints_mgi_only,"imprints_from_literature_mgi_only.csv",row.names=FALSE)



### testing if some genes are dropped due to cutoffs or not found to be imprinted ####
## list of known imprinted genes from Blake, A. et al., Nucleic Acids 2010 (doi:10.1093/nar/gkp867)
## https://www.mousebook.org/
## not provided in our upload
impr_lit <- read.xls("mousebook imprinted genes.xlsx")
### 76 imprinted genes from Inoue et al., Nature 2017 (https://www.nature.com/articles/nature23262)
### file is not provided in this upload
inoue <- read.xls("Inoue_76.xlsx")

## file from other scripts of the analysis
Florian_genes <- read.csv("d_disjointexons_counts.csv",
                          stringsAsFactors = F)

## files retrieved from featurCounts of mapped BAM files
all_counts <- data.frame(names=as.character(read.table("B6_x_cast_1xE3_5_bla_1.txt")[1:106519,1]),
                         mgi=NA,
                         B6xcast_1=read.table("B6_x_cast_1xE3_5_bla_1.txt")[1:106519,2],
                         B6xcast_2=read.table("B6_x_cast_1xE3_5_bla_2.txt")[1:106519,2],
                         B6xcast_3=read.table("B6_x_cast_1xE3_5_bla_3.txt")[1:106519,2],
                         castxB6_1=read.table("cast_x_B6_1xE3_5_bla_1.txt")[1:106519,2],
                         castxB6_2=read.table("cast_x_B6_1xE3_5_bla_2.txt")[1:106519,2],
                         castxB6_3=read.table("cast_x_B6_1xE3_5_bla_3.txt")[1:106519,2],
                         castxB6_5=read.table("cast_x_B6_1xE3_5_bla_5.txt")[1:106519,2],
                         castxB6_7=read.table("cast_x_B6_1xE3_5_bla_7.txt")[1:106519,2])

dim(all_counts[which(!all_counts$names%in%ncbiRef$V2),])

for(d in 1:dim(all_counts)[1]){
  if(d%%10000==0){
    print(d)
  }
  all_counts$mgi[d] <- ncbiRef[which(ncbiRef$V2==all_counts$names[d]),2]
}

robustly_expressed <- Florian_genes[which(apply(Florian_genes[,2:9],1, function(x) sum(x >= 12) >= 4)),1]

length(unique(all_counts[which(all_counts$mgi%in%inoue$Gene.symbol),2]))
#76

dim(merged_counts[which(row.names(merged_counts)%in%inoue$Gene.symbol),])
#48

dim(counts_mx[which(row.names(counts_mx)%in%inoue$Gene.symbol),])
#36

length(unique(ncbiRef[which(ncbiRef$V13%in%inoue$Gene.symbol),2]))
#76

dim(Florian_genes[which(Florian_genes$rn%in%inoue$Gene.symbol),])
#36

length(robustly_expressed[which(robustly_expressed%in%inoue$Gene.symbol)])
#29

## file used for Allelome.Pro analysis
## provided in this upload
snp_genes <- read.table("anno_bigger100_nonoverlap.bed")
snp_genes <- unique(snp_genes$V4)
sum(snp_genes%in%inoue$Gene.symbol)
#76

### control gene groups in TC from Lackner et al 2021 ####
## HC publicly available imprinted genes
## provided in this upload
HC_genes <- as.character(read.csv("HC_repositories_imprints.csv",header = F)[,1])

## table containing gene groups based on new and known imprints as well as skewed or imprinted genes
## provided in this upload
final_table_laura <- openxlsx::read.xlsx("All_detected_genes_FINAL.xlsx",sheet = 1,startRow=2)


Bix_LS <- final_table_laura[which(final_table_laura$Is.Bix.LS==1),1]
nBix_LS <- final_table_laura[which(final_table_laura$Is.nBix.LS==1),1]

Bix_HS <- final_table_laura[which(final_table_laura$Is.Bix.HS==1),1]
nBix_HS <- final_table_laura[which(final_table_laura$Is.nBix.HS==1),1]

published_imprint <- final_table_laura[which(final_table_laura$Published.imprint==1),1]

## read in TC data from Lackner et al 2021
## provided in this upload
table_boxplots_prep <- readRDS("TC_table.rds")



ggplot_tidy_frame_24_nonabs <- data.frame(names=table_boxplots_prep[which(table_boxplots_prep$max>1),24],
                                          log2FCs=table_boxplots_prep[which(table_boxplots_prep$max>1),17],
                                          group="all expressed genes")

ggplot_tidy_frame_24_nonabs <- rbind(ggplot_tidy_frame_24_nonabs,
                                     data.frame(names=table_boxplots_prep[which(table_boxplots_prep$mgi%in%HC_genes),24],
                                                log2FCs=table_boxplots_prep[which(table_boxplots_prep$mgi%in%HC_genes),17],
                                                group="HC genes"))

ggplot_tidy_frame_24_nonabs <- rbind(ggplot_tidy_frame_24_nonabs,
                                     data.frame(names=table_boxplots_prep[which(table_boxplots_prep$mgi%in%nBix_HS),24],
                                                log2FCs=table_boxplots_prep[which(table_boxplots_prep$mgi%in%nBix_HS),17],
                                                group="nBix"))

ggplot_tidy_frame_24_nonabs <- rbind(ggplot_tidy_frame_24_nonabs,
                                     data.frame(names=table_boxplots_prep[which(table_boxplots_prep$mgi%in%nBix_LS),24],
                                                log2FCs=table_boxplots_prep[which(table_boxplots_prep$mgi%in%nBix_LS),17],
                                                group="nBsx"))

ggplot_tidy_frame_24_nonabs <- rbind(ggplot_tidy_frame_24_nonabs,
                                     data.frame(names=table_boxplots_prep[which(table_boxplots_prep$mgi%in%Bix_HS),24],
                                                log2FCs=table_boxplots_prep[which(table_boxplots_prep$mgi%in%Bix_HS),17],
                                                group="Bix"))

ggplot_tidy_frame_24_nonabs <- rbind(ggplot_tidy_frame_24_nonabs,
                                     data.frame(names=table_boxplots_prep[which(table_boxplots_prep$mgi%in%Bix_LS),24],
                                                log2FCs=table_boxplots_prep[which(table_boxplots_prep$mgi%in%Bix_LS),17],
                                                group="Bsx"))

ggplot_tidy_frame_24_nonabs <- rbind(ggplot_tidy_frame_24_nonabs,
                                     data.frame(names=table_boxplots_prep[which(table_boxplots_prep$mgi%in%published_imprint),24],
                                                log2FCs=table_boxplots_prep[which(table_boxplots_prep$mgi%in%published_imprint),17],
                                                group="published imprint"))

ggplot_tidy_frame_24_nonabs <- ggplot_tidy_frame_24_nonabs[which(ggplot_tidy_frame_24_nonabs$names%in%robustly_expressed),]

ggplot(data=ggplot_tidy_frame_24_nonabs[which(ggplot_tidy_frame_24_nonabs$group%in%
                                                c("all expressed genes",
                                                  "published imprint","HC genes",
                                                  "nBsx","nBix",
                                                  "Bsx","Bix")),],
       aes(x=factor(group,levels=c("all expressed genes",
                                   "published imprint","HC genes",
                                   "nBsx","nBix",
                                   "Bsx","Bix")),
       y=log2FCs)) +
  geom_violin(draw_quantiles = 0.5) +
  ylab("log2FC 24 hours vs 2i") +
  xlab("groups") +
  ggtitle("Comparison of log2FCs for different groups of genes\n") +
  theme(plot.title=element_text(hjust=0.5)) + coord_flip() 


#write.csv(ggplot_tidy_frame_24_nonabs, "violin_data_050321.csv", row.names = F, quote = F)