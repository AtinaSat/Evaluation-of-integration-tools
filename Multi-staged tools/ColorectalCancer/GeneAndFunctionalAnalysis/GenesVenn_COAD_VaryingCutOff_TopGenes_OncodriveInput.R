## Varying cutoff Gene-level comparison of results from copy number and gene expression integration tools for colorectal cancer dataset
# with cBioPortal cancer gene list

# In this script, we read the results from copy number and gene integration tools (CNAmet, iGC, plrs and Oncodrive-CIS)
# and select the top 3000 genes based on fdr values and perform gene-level comparison with cancer gene list from cBioPortal
# cBioPortal cancer genes is predownloaded and saved as cBioPortalCancerGeneList.txt

## Cut-offs are at top 1000 (500 amp, 500 del), top 2000 (1000 amp, 1000 del), and top 3000 (1500 amp, 1500 del)


rm(list=ls())

## CNAmet does not perform fdr correction if only two omics data are integrated. 
# hence, perform fdr correction and select top 1500 amplification driven and top 1500 deletion driven genes based on fdr

# amplification driven genes
CNAmetAmp <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_COAD_GainDriven_Genes_OncodriveInput.tsv",
                        header = T, sep = "\t") # 12358
CNAmetAmp <- CNAmetAmp[,c(3,4,10)]
CNAmetAmp <- na.omit(CNAmetAmp) # 10262
CNAmetAmp$p.adj.fdr <- p.adjust(CNAmetAmp$CWPvalue, method = "fdr")
CNAmetAmp <- CNAmetAmp[CNAmetAmp$p.adj.fdr<0.05,] # 6471
CNAmetAmp <- CNAmetAmp[order(CNAmetAmp$p.adj.fdr),]
CNAmetAmpSelected1000 <- CNAmetAmp[1:500,]
CNAmetAmpSelected2000 <- CNAmetAmp[1:1000,]
CNAmetAmpSelected3000 <- CNAmetAmp[1:1500,]
CNAmetAmpSelected <- CNAmetAmp

# deletion driven genes
CNAmetDel <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_COAD_LossDriven_Genes_OncodriveInput.tsv",
                        header = T, sep = "\t") # 12358
CNAmetDel <- CNAmetDel[,c(3,4,10)]
CNAmetDel <- na.omit(CNAmetDel) # 9614
CNAmetDel$p.adj.fdr <- p.adjust(CNAmetDel$CWPvalue, method = "fdr")
CNAmetDel <- CNAmetDel[CNAmetDel$p.adj.fdr<0.05,] # 5818
CNAmetDel <- CNAmetDel[order(CNAmetDel$p.adj.fdr),]
CNAmetDelSelected1000 <- CNAmetDel[1:500,]
CNAmetDelSelected2000 <- CNAmetDel[1:1000,]
CNAmetDelSelected3000 <- CNAmetDel[1:1500,]
CNAmetDelSelected <- CNAmetDel

# Both genes
CNAmetBoth1000 <- intersect(as.character(CNAmetAmpSelected1000$Gene), as.character(CNAmetDelSelected1000$Gene)) # 69
CNAmetBoth2000 <- intersect(as.character(CNAmetAmpSelected2000$Gene), as.character(CNAmetDelSelected2000$Gene)) # 233
CNAmetBoth3000 <- intersect(as.character(CNAmetAmpSelected3000$Gene), as.character(CNAmetDelSelected3000$Gene)) # 461
CNAmetBoth <- intersect(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene)) # 4205

# Total genes
CNAmetTotal1000 <- unique(c(as.character(CNAmetAmpSelected1000$Gene), as.character(CNAmetDelSelected1000$Gene))) # 931
CNAmetTotal2000 <- unique(c(as.character(CNAmetAmpSelected2000$Gene), as.character(CNAmetDelSelected2000$Gene))) # 1767
CNAmetTotal3000 <- unique(c(as.character(CNAmetAmpSelected3000$Gene), as.character(CNAmetDelSelected3000$Gene))) # 2539
CNAmetTotal <- unique(c(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene))) # 8084

# save genes
write.table(CNAmetTotal1000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_COAD_1000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(CNAmetTotal2000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_COAD_2000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(CNAmetTotal3000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_COAD_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

write.table(CNAmetTotal, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_COAD_AllGenes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################

## iGC
# read the results and select top 1500 amplification driven and 1500 deletion driven genes based on fdr

# amplification driven genes
iGCAmp <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_COAD_GainDriven_genes_OncodriveInput.tsv",
                     header = T, sep = "\t") # 4990
iGCAmp <- iGCAmp[,c(1,3)]
iGCAmp <- iGCAmp[iGCAmp$fdr<0.05,] # 2870


# deletion driven genes
iGCDel <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_COAD_LossDriven_genes_OncodriveInput.tsv",
                     header = T, sep = "\t") # 3499
iGCDel <- iGCDel[,c(1,3)]
iGCDel <- iGCDel[iGCDel$fdr<0.05,] # 2231

# Both genes
iGCBoth <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_COAD_BothDriven_genes_OncodriveInput.tsv",
                      header = T, sep = "\t") # 236
iGCBothAmp <- iGCBoth[(iGCBoth$gain_fdr<0.05),] # 174
iGCBothAmp <- iGCBothAmp[,c(1,3)]
colnames(iGCBothAmp)[2] <- "fdr"

iGCBothDel <- iGCBoth[(iGCBoth$loss_fdr<0.05),] # 177
iGCBothDel <- iGCBothDel[,c(1,5)]
colnames(iGCBothDel)[2] <- "fdr"

iGCAmp <- rbind(iGCAmp, iGCBothAmp) # 3044
iGCAmp <- iGCAmp[order(iGCAmp$fdr),]
iGCAmpSelected1000 <- iGCAmp[1:500,]
iGCAmpSelected2000 <- iGCAmp[1:1000,]
iGCAmpSelected3000 <- iGCAmp[1:1500,]
iGCAmpSelected <- iGCAmp

iGCDel <- rbind(iGCDel, iGCBothDel)
iGCDel <- iGCDel[order(iGCDel$fdr),]
iGCDelSelected1000 <- iGCDel[1:500,]
iGCDelSelected2000 <- iGCDel[1:1000,]
iGCDelSelected3000 <- iGCDel[1:1500,]
iGCDelSelected <- iGCDel


# Both and Total genes
iGCBoth1000 <- intersect(as.character(iGCAmpSelected1000$GENE), as.character(iGCDelSelected1000$GENE)) # 70
iGCBoth2000 <- intersect(as.character(iGCAmpSelected2000$GENE), as.character(iGCDelSelected2000$GENE)) # 111
iGCBoth3000 <- intersect(as.character(iGCAmpSelected3000$GENE), as.character(iGCDelSelected3000$GENE)) # 121
iGCBoth <- intersect(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE)) # 167

iGCTotal1000 <- unique(c(as.character(iGCAmpSelected1000$GENE), as.character(iGCDelSelected1000$GENE))) # 930
iGCTotal2000 <- unique(c(as.character(iGCAmpSelected2000$GENE), as.character(iGCDelSelected2000$GENE))) # 1889
iGCTotal3000 <- unique(c(as.character(iGCAmpSelected3000$GENE), as.character(iGCDelSelected3000$GENE))) # 2879
iGCTotal <- unique(c(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE))) # 5285

# save genes
write.table(iGCTotal1000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_COAD_1000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(iGCTotal2000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_COAD_2000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(iGCTotal3000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_COAD_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(iGCTotal, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_COAD_AllGenes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## plrs
# plrs does not output amplfied and deletion driven genes. Hence selected top 3000 genes from results based on fdr

plrsGenes <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_COAD_Results_OncodriveInput.tsv",
                        header = T, sep = "\t") # 12358
plrsGenes <- plrsGenes[plrsGenes$BH.adj.pval<0.05,]  # 8039
plrsGenes <- plrsGenes[order(plrsGenes$BH.adj.pval),]
plrsSelected1000 <- plrsGenes[1:1000,]
plrsSelected2000 <- plrsGenes[1:2000,]
plrsSelected3000 <- plrsGenes[1:3000,]
plrsSelected <- plrsGenes

plrsTotal1000 <- as.character(plrsSelected1000$Gene)
plrsTotal2000 <- as.character(plrsSelected2000$Gene)
plrsTotal3000 <- as.character(plrsSelected3000$Gene)
plrsTotal <- as.character(plrsSelected$Gene)

# save genes
write.table(plrsTotal1000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_COAD_1000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(plrsTotal2000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_COAD_2000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(plrsTotal3000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_COAD_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(plrsTotal, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_COAD_AllGenes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## Oncodrive-CIS
# did not provide p-values. Hence we estimate p-values using z-score and estimate fdr
# select top 1500 amplification driven and 1500 deletion driven gene based on fdr

# amplification driven genes
OncoAmp <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/OncodriveCIS/CompCOAD.Output/OncoCNA.AMP",
                      header = T, sep = "\t") # 12358
OncoAmp$pval <- pnorm(-abs(OncoAmp$Z.COMBINED))
OncoAmp$p.adj.fdr <- p.adjust(OncoAmp$pval, method = "fdr")
OncoAmp <- OncoAmp[OncoAmp$Z.COMBINED>0,]
OncoAmp <- OncoAmp[OncoAmp$p.adj.fdr<0.05,] # 1518
OncoAmp <- OncoAmp[order(OncoAmp$p.adj.fdr),]
OncoAmpSelected1000 <- OncoAmp[1:500,]
OncoAmpSelected2000 <- OncoAmp[1:1000,]
OncoAmpSelected3000 <- OncoAmp[1:1500,]
OncoAmpSelected <- OncoAmp

# deletion driven genes
OncoDel <- read.delim("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/OncodriveCIS/CompCOAD.Output/OncoCNA.DEL",
                      header = T, sep = "\t") # 11891
head(OncoDel)
tail(OncoDel)
OncoDel$pval <- pnorm(-abs(OncoDel$Z.COMBINED))
OncoDel$p.adj.fdr <- p.adjust(OncoDel$pval, method = "fdr")
OncoDel <- OncoDel[OncoDel$Z.COMBINED<0,]
OncoDel <- OncoDel[OncoDel$p.adj.fdr<0.05,] # 1771
OncoDel <- OncoDel[order(OncoDel$p.adj.fdr),]
OncoDelSelected1000 <- OncoDel[1:500,]
OncoDelSelected2000 <- OncoDel[1:1000,]
OncoDelSelected3000 <- OncoDel[1:1500,]
OncoDelSelected <- OncoDel

# Both genes
OncoBoth1000 <- intersect(as.character(OncoAmpSelected1000$ENS.ID), as.character(OncoDelSelected1000$ENS.ID)) # 30
OncoBoth2000 <- intersect(as.character(OncoAmpSelected2000$ENS.ID), as.character(OncoDelSelected2000$ENS.ID)) # 115
OncoBoth3000 <- intersect(as.character(OncoAmpSelected3000$ENS.ID), as.character(OncoDelSelected3000$ENS.ID)) # 267
OncoBoth <- intersect(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID)) # 332

# Total genes
OncoTotal1000 <- unique(c(as.character(OncoAmpSelected1000$ENS.ID), as.character(OncoDelSelected1000$ENS.ID))) # 970
OncoTotal2000 <- unique(c(as.character(OncoAmpSelected2000$ENS.ID), as.character(OncoDelSelected2000$ENS.ID))) # 1885
OncoTotal3000 <- unique(c(as.character(OncoAmpSelected3000$ENS.ID), as.character(OncoDelSelected3000$ENS.ID))) # 2733
OncoTotal <- unique(c(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID))) # 2957

# save genes
write.table(OncoTotal1000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/OncodriveCIS/Oncodrive_UsingOncodriveCISInput/Oncodrive_COAD_1000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(OncoTotal2000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/OncodriveCIS/Oncodrive_UsingOncodriveCISInput/Oncodrive_COAD_2000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(OncoTotal3000, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/OncodriveCIS/Oncodrive_UsingOncodriveCISInput/Oncodrive_COAD_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")
write.table(OncoTotal, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/OncodriveCIS/Oncodrive_UsingOncodriveCISInput/Oncodrive_COAD_AllGenes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## cBioPortal Cancer Gene List and Venn
# copied the list of genes from http://www.cbioportal.org/cancer_gene_list.jsp onto a text file

library(readr)
# reading the text file
mystring <- read_file("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/Venns/cBioPortalCancerGeneList.txt")
cancergenes <- strsplit(mystring, "\n", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]

CommonGenesVarying <- data.frame(matrix(ncol = 5, nrow = 4))
colnames(CommonGenesVarying) <- c("Tool", "Top 1000", "Top 2000", "Top 3000", "All")
CommonGenesVarying$Tool <- c("CNAmet", "iGC", "PLRS", "OncodriveCIS")

# finding common genes between tools and cancer gene list
CommonGenesVarying[1,2] <- length(cancerCNAmet1000 <- intersect(cancergenes, CNAmetTotal1000))
CommonGenesVarying[2,2] <- length(canceriGC1000 <- intersect(cancergenes, iGCTotal1000))
CommonGenesVarying[3,2] <- length(cancerplrs1000 <- intersect(cancergenes, plrsTotal1000))
CommonGenesVarying[4,2] <- length(cancerOnco1000 <- intersect(cancergenes, OncoTotal1000))

# finding common genes between tools and cancer gene list
CommonGenesVarying[1,3] <- length(cancerCNAmet2000 <- intersect(cancergenes, CNAmetTotal2000))
CommonGenesVarying[2,3] <- length(canceriGC2000 <- intersect(cancergenes, iGCTotal2000))
CommonGenesVarying[3,3] <- length(cancerplrs2000 <- intersect(cancergenes, plrsTotal2000))
CommonGenesVarying[4,3] <- length(cancerOnco2000 <- intersect(cancergenes, OncoTotal2000))

# finding common genes between tools and cancer gene list
CommonGenesVarying[1,4] <- length(cancerCNAmet3000 <- intersect(cancergenes, CNAmetTotal3000))
CommonGenesVarying[2,4] <- length(canceriGC3000 <- intersect(cancergenes, iGCTotal3000))
CommonGenesVarying[3,4] <- length(cancerplrs3000 <- intersect(cancergenes, plrsTotal3000))
CommonGenesVarying[4,4] <- length(cancerOnco3000 <- intersect(cancergenes, OncoTotal3000))

# finding common genes between tools and cancer gene list
CommonGenesVarying[1,5] <- length(cancerCNAmet <- intersect(cancergenes, CNAmetTotal))
CommonGenesVarying[2,5] <- length(canceriGC <- intersect(cancergenes, iGCTotal))
CommonGenesVarying[3,5] <- length(cancerplrs <- intersect(cancergenes, plrsTotal))
CommonGenesVarying[4,5] <- length(cancerOnco <- intersect(cancergenes, OncoTotal))


setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/Venns_OncodriveInput/")
write.table(CommonGenesVarying, file = "COAD_VaryingCutOff_TopGenes_OncodriveInput_Results.txt", sep = "\t", row.names = F, quote = F)



# drawing venn diagram
library(venn)
venn(list(CNAmet = CNAmetTotal,  iGC = iGCTotal,  PLRS = plrsTotal,
          OncodriveCIS = OncoTotal,
          CCGL = cancergenes),
     zcolor = "style", cexil = 1.5, cexsn = 1.5, opacity = 0.3,size = 20)


## Save the venn plot
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/Venns_OncodriveInput//")

tiff("COAD_OncodriveInput_GenesVenn.tiff", width = 952, height = 952, res = 300)  
venn(list(CNAmet = CNAmetTotal,  iGC = iGCTotal,  PLRS = plrsTotal,
          OncodriveCIS = OncoTotal,
          CCGL = cancergenes),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
