## Gene-level comparison of results from copy number and gene expression integration tools for Mesothelioma dataset using 
# Oncodrive CIS Input
# with cBioPortal cancer gene list 
# Date: 26/07/2019

# In this script, we read the results from copy number and gene integration tools (CNAmet, iGC, plrs and Oncodrive-CIS)
# and select the top 3000 genes based on fdr values 

rm(list=ls())

## CNAmet does not perform fdr correction if only two omics data are integrated. 
# hence, perform fdr correction and select top 1500 amplification driven and top 1500 deletion driven genes based on fdr

# amplification driven genes
CNAmetAmp <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_MESO_GainDriven_Genes_OncodriveInput.tsv",
                        header = T, sep = "\t") # 8372
CNAmetAmp <- CNAmetAmp[,c(3,4,10)]
CNAmetAmp <- na.omit(CNAmetAmp)
CNAmetAmp$p.adj.fdr <- p.adjust(CNAmetAmp$CWPvalue, method = "fdr")
CNAmetAmp <- CNAmetAmp[CNAmetAmp$p.adj.fdr<0.05,] # 719 genes
CNAmetAmp <- CNAmetAmp[order(CNAmetAmp$p.adj.fdr),]
CNAmetAmpSelected <- CNAmetAmp

# deletion driven genes
CNAmetDel <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_MESO_LossDriven_Genes_OncodriveInput.tsv",
                        header = T, sep = "\t") # 8372
CNAmetDel <- CNAmetDel[,c(3,4,10)]
CNAmetDel <- na.omit(CNAmetDel)
CNAmetDel$p.adj.fdr <- p.adjust(CNAmetDel$CWPvalue, method = "fdr")
CNAmetDel <- CNAmetDel[CNAmetDel$p.adj.fdr<0.05,] # 1815
CNAmetDel <- CNAmetDel[order(CNAmetDel$p.adj.fdr),]
CNAmetDelSelected <- CNAmetDel

# Both genes
CNAmetBoth <- intersect(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene)) # 262

# Total genes
CNAmetTotal <- unique(c(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene))) # 2272

# save genes
write.table(CNAmetTotal, file = "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_MESO_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## iGC
# read the results and select top 1500 amplification driven and 1500 deletion driven genes based on fdr

# amplification driven genes
iGCAmp <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_MESO_GainDriven_genes_OncodriveInput.tsv",
                     header = T, sep = "\t") # 1819
iGCAmp <- iGCAmp[,c(1,3)]
iGCAmp <- iGCAmp[iGCAmp$fdr<0.05,] # 231

# deletion driven genes
iGCDel <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_MESO_LossDriven_genes_OncodriveInput.tsv",
                     header = T, sep = "\t") # 2630
iGCDel <- iGCDel[,c(1,3)]
iGCDel <- iGCDel[iGCDel$fdr<0.05,] # 924

# Both genes
iGCBoth <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_MESO_BothDriven_genes_OncodriveInput.tsv",
                      header = T, sep = "\t") # 0
# There are no genes in iGCBoth and hence we don't run the following till line 73
# iGCBothAmp <- iGCBoth[(iGCBoth$gain_fdr<0.05),]
# iGCBothAmp <- iGCBothAmp[,c(1,3)]
# colnames(iGCBothAmp)[2] <- "fdr"
# 
# iGCBothDel <- iGCBoth[(iGCBoth$loss_fdr<0.05),]
# iGCBothDel <- iGCBothDel[,c(1,5)]
# colnames(iGCBothDel)[2] <- "fdr"
# 
# iGCAmp <- rbind(iGCAmp, iGCBothAmp)
# iGCAmpSelected <- iGCAmp # less than 1500
# iGCDel <- rbind(iGCDel, iGCBothDel)
# iGCDel <- iGCDel[order(iGCDel$fdr),]

iGCDel <- iGCDel[order(iGCDel$fdr),]
iGCDelSelected <- iGCDel # less than 1500
iGCAmpSelected <- iGCAmp # less than 1500

# Both genes and Total Genes
iGCBoth <- intersect(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE)) # 0
iGCTotal <- unique(c(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE))) # 1155

# save genes
write.table(iGCTotal, file = "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_MESO_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## plrs
# plrs does not output amplfied and deletion driven genes. Hence selected top 3000 genes from results based on fdr

# total genes
plrsGenes <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_MESO_Results_OncodriveInput.tsv",
                        header = T, sep = "\t") # 8372
plrsGenes<- plrsGenes[plrsGenes$BH.adj.pval<0.05,] # 2472
plrsGenes <- plrsGenes[order(plrsGenes$BH.adj.pval),]
plrsSelected <- plrsGenes # less than 3000

plrsTotal <- as.character(plrsSelected$Gene)

# save genes
write.table(plrsTotal, file = "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_MESO_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## Oncodrive-CIS
# did not provide p-values. Hence we estimate p-values using z-score and estimate fdr
# select top 1500 amplification driven and 1500 deletion driven gene based on fdr

# amplification driven genes
OncoAmp <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/CompMESO.Output/OncoCNA.AMP",
                      header = T, sep = "\t") # 7861
OncoAmp$pval <- pnorm(-abs(OncoAmp$Z.COMBINED))
OncoAmp$p.adj.fdr <- p.adjust(OncoAmp$pval, method = "fdr")
OncoAmp <- OncoAmp[OncoAmp$Z.COMBINED>0,]
OncoAmpSelected <- OncoAmp[OncoAmp$p.adj.fdr<0.05,] # less than 1500.. so using all # 144 genes

# deletion driven genes
OncoDel <- read.delim("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/CompMESO.Output/OncoCNA.DEL",
                      header = T, sep = "\t") # 7976
head(OncoDel)
tail(OncoDel)
OncoDel$pval <- pnorm(-abs(OncoDel$Z.COMBINED))
OncoDel$p.adj.fdr <- p.adjust(OncoDel$pval, method = "fdr")
OncoDel <- OncoDel[OncoDel$Z.COMBINED<0,]
OncoDelSelected <- OncoDel[OncoDel$p.adj.fdr<0.05,] # less than 1500 .. so using all # 532 genes

# Both genes
OncoBoth <- intersect(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID)) # 13 genes

# Total genes
OncoTotal <- unique(c(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID))) # 663

# save genes
write.table(OncoTotal, file = "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/OncodriveCIS/Oncodrive_UsingOncodriveCISInput/Oncodrive_MESO_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")


####################################################################################################################
## cBioPortal Cancer Gene List and Venn
# copied the list of genes from http://www.cbioportal.org/cancer_gene_list.jsp onto a text file

library(readr)
# reading the text file
mystring <- read_file("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/Venns/cBioPortalCancerGeneList.txt")
cancergenes <- strsplit(mystring, "\n", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]

# finding common genes between tools and cancer gene list
cancerCNAmet <- intersect(cancergenes, CNAmetTotal)
canceriGC <- intersect(cancergenes, iGCTotal)
cancerplrs <- intersect(cancergenes, plrsTotal)
cancerOnco <- intersect(cancergenes, OncoTotal)

# drawing venn diagram
library(venn)
venn(list(CNAmet = CNAmetTotal,  iGC = iGCTotal,  PLRS = plrsTotal,
          OncodriveCIS = OncoTotal,
          CCGL = cancergenes),
     zcolor = "style", cexil = 1.5, cexsn = 1.5, opacity = 0.3,size = 20)

## Save the venn plot
setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/Venns_OncodriveInput/")
tiff("MESO_OncodriveInput_GenesVenn.tiff", width = 952, height = 952, res=300)  
venn(list(CNAmet = CNAmetTotal,  iGC = iGCTotal,  PLRS = plrsTotal,
          OncodriveCIS = OncoTotal,
          CCGL = cancergenes),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()

