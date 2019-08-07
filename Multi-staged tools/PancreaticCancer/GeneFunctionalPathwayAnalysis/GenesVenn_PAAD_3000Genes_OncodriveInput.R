## Gene-level comparison of results from copy number and gene expression integration tools for Pancreatic cancer dataset
# with cBioPortal cancer gene list

# In this script, we read the results from copy number and gene integration tools (CNAmet, iGC, plrs and Oncodrive-CIS)
# and select the top 3000 genes based on fdr values and perform gene-level comparison with cancer gene list from cBioPortal
# cBioPortal cancer genes is predownloaded and saved as cBioPortalCancerGeneList.txt

rm(list=ls())

## CNAmet does not perform fdr correction if only two omics data are integrated. 
# hence, perform fdr correction and select top 1500 amplification driven and top 1500 deletion driven genes based on fdr

# amplification driven genes
CNAmetAmp <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_PAAD_GainDriven_Genes_OncodriveInput.tsv",
                        header = T, sep = "\t") # 12977
CNAmetAmp <- CNAmetAmp[,c(3,4,10)]
CNAmetAmp <- na.omit(CNAmetAmp) # 8695
CNAmetAmp$p.adj.fdr <- p.adjust(CNAmetAmp$CWPvalue, method = "fdr")
CNAmetAmp <- CNAmetAmp[CNAmetAmp$p.adj.fdr<0.05,] # 3927
CNAmetAmp <- CNAmetAmp[order(CNAmetAmp$p.adj.fdr),]
CNAmetAmpSelected <- CNAmetAmp

# deletion driven genes
CNAmetDel <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_PAAD_LossDriven_Genes_OncodriveInput.tsv",
                        header = T, sep = "\t") # 12977
CNAmetDel <- CNAmetDel[,c(3,4,10)]
CNAmetDel <- na.omit(CNAmetDel) # 11557
CNAmetDel$p.adj.fdr <- p.adjust(CNAmetDel$CWPvalue, method = "fdr")
CNAmetDel <- CNAmetDel[CNAmetDel$p.adj.fdr<0.05,] # 6874
CNAmetDel <- CNAmetDel[order(CNAmetDel$p.adj.fdr),]
CNAmetDelSelected <- CNAmetDel

# Both genes
CNAmetBoth <- intersect(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene))

# Total genes
CNAmetTotal <- unique(c(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene)))

# save genes
write.table(CNAmetTotal, file = "/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_PAAD_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## iGC
# read the results and select top 1500 amplification driven and 1500 deletion driven genes based on fdr

# amplification driven genes
iGCAmp <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_PAAD_GainDriven_genes_OncodriveInput.tsv",
                     header = T, sep = "\t") # 2199
iGCAmp <- iGCAmp[,c(1,3)]
iGCAmp <- iGCAmp[iGCAmp$fdr<0.05,] # 988

# deletion driven genes
iGCDel <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_PAAD_LossDriven_genes_OncodriveInput.tsv",
                     header = T, sep = "\t") # 4027
iGCDel <- iGCDel[,c(1,3)]
iGCDel <- iGCDel[iGCDel$fdr<0.05,] # 2467

# Both genes
iGCBoth <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_PAAD_BothDriven_genes_OncodriveInput.tsv",
                      header = T, sep = "\t") # 16
iGCBothAmp <- iGCBoth[(iGCBoth$gain_fdr<0.05),] # 12
iGCBothAmp <- iGCBothAmp[,c(1,3)]
colnames(iGCBothAmp)[2] <- "fdr"

iGCBothDel <- iGCBoth[(iGCBoth$loss_fdr<0.05),] # 9
iGCBothDel <- iGCBothDel[,c(1,5)]
colnames(iGCBothDel)[2] <- "fdr"

iGCAmp <- rbind(iGCAmp, iGCBothAmp) # 1000
iGCAmpSelected <- iGCAmp # less than 1500
iGCDel <- rbind(iGCDel, iGCBothDel) # 2476
iGCDel <- iGCDel[order(iGCDel$fdr),]
iGCDelSelected <- iGCDel

# Both genes and Total Genes
iGCBoth <- intersect(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE))
iGCTotal <- unique(c(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE)))

# save genes
write.table(iGCTotal, file = "/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_PAAD_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## plrs
# plrs does not output amplfied and deletion driven genes. Hence selected top 3000 genes from results based on fdr

# total genes
plrsGenes <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_PAAD_Results_OncodriveInput.tsv",
                        header = T, sep = "\t") # 12977
plrsGenes<- plrsGenes[plrsGenes$BH.adj.pval<0.05,]  # 8369
plrsGenes <- plrsGenes[order(plrsGenes$BH.adj.pval),]
plrsSelected <- plrsGenes # 8369
plrsTotal <- as.character(plrsSelected$Gene)

# save genes
write.table(plrsTotal, file = "/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_PAAD_3000Genes_OncodriveInput.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## Oncodrive-CIS
# did not provide p-values. Hence we estimate p-values using z-score and estimate fdr
# select top 1500 amplification driven and 1500 deletion driven gene based on fdr

# amplification driven genes
OncoAmp <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/OncodriveCIS/CompPAAD.Output/OncoCNA.AMP",
                      header = T, sep = "\t") # 12977
OncoAmp$pval <- pnorm(-abs(OncoAmp$Z.COMBINED))
OncoAmp$p.adj.fdr <- p.adjust(OncoAmp$pval, method = "fdr")
OncoAmp <- OncoAmp[OncoAmp$Z.COMBINED>0,]
OncoAmpSelected <- OncoAmp[OncoAmp$p.adj.fdr<0.05,] # less than 1500.. so using all # 1279

# deletion driven genes
OncoDel <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/OncodriveCIS/CompPAAD.Output/OncoCNA.DEL",
                      header = T, sep = "\t") # 12977
head(OncoDel)
tail(OncoDel)
OncoDel$pval <- pnorm(-abs(OncoDel$Z.COMBINED))
OncoDel$p.adj.fdr <- p.adjust(OncoDel$pval, method = "fdr")
OncoDel <- OncoDel[OncoDel$Z.COMBINED<0,]
OncoDelSelected <- OncoDel[OncoDel$p.adj.fdr<0.05,] # less than 1500 .. so using all # 1048

# Both genes
OncoBoth <- intersect(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID))

# Total genes
OncoTotal <- unique(c(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID)))

# save genes
write.table(OncoTotal, file = "/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/OncodriveCIS/Oncodrive_PAAD_3000Genes.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## cBioPortal Cancer Gene List and Venn
# copied the list of genes from http://www.cbioportal.org/cancer_gene_list.jsp onto a text file

library(readr)
# reading the text file
mystring <- read_file("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns/cBioPortalCancerGeneList.txt")
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
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns_OncodriveInput/")
tiff("PAAD_OncodriveInput_GenesVenn.tiff", width = 952, height = 952, res=300)  
venn(list(CNAmet = CNAmetTotal,  iGC = iGCTotal,  PLRS = plrsTotal,
          OncodriveCIS = OncoTotal,
          CCGL = cancergenes),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()

