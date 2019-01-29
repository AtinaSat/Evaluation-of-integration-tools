## Gene-level comparison of results from copy number and gene expression integration tools for Pancreatic cancer dataset
# with cBioPortal cancer gene list

# In this script, we read the results from copy number and gene integration tools (CNAmet, iGC, plrs and Oncodrive-CIS)
# and select the top 3000 genes based on fdr values and perform gene-level comparison with cancer gene list from cBioPortal
# cBioPortal cancer genes is predownloaded and saved as cBioPortalCancerGeneList.txt

rm(list=ls())

## CNAmet does not perform fdr correction if only two omics data are integrated. 
# hence, perform fdr correction and select top 1500 amplification driven and top 1500 deletion driven genes based on fdr

# amplification driven genes
CNAmetAmp <- read.delim("PancreaticCancer/CNAmet/CNAmet_PAAD_GainDriven_Genes.tsv", header = T, sep = "\t")
CNAmetAmp <- CNAmetAmp[,c(3,4,10)]
CNAmetAmp <- na.omit(CNAmetAmp)
CNAmetAmp$p.adj.fdr <- p.adjust(CNAmetAmp$CWPvalue, method = "fdr")
CNAmetAmp <- CNAmetAmp[CNAmetAmp$p.adj.fdr<0.05,]
CNAmetAmp <- CNAmetAmp[order(CNAmetAmp$p.adj.fdr),]
CNAmetAmpSelected <- CNAmetAmp[1:1500,]

# deletion driven genes
CNAmetDel <- read.delim("PancreaticCancer/CNAmet/CNAmet_PAAD_LossDriven_Genes.tsv",header = T, sep = "\t")
CNAmetDel <- CNAmetDel[,c(3,4,10)]
CNAmetDel <- na.omit(CNAmetDel)
CNAmetDel$p.adj.fdr <- p.adjust(CNAmetDel$CWPvalue, method = "fdr")
CNAmetDel <- CNAmetDel[CNAmetDel$p.adj.fdr<0.05,]
CNAmetDel <- CNAmetDel[order(CNAmetDel$p.adj.fdr),]
CNAmetDelSelected <- CNAmetDel[1:1500,]

# Both genes
CNAmetBoth <- intersect(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene))

# Total genes
CNAmetTotal <- unique(c(as.character(CNAmetAmpSelected$Gene), as.character(CNAmetDelSelected$Gene)))

# save genes
# write.table(CNAmetTotal, file = "CNAmet/CNAmet_PAAD_3000Genes.txt", row.names = F, quote = F, sep = "\t")

####################################################################################################################
## iGC
# read the results and select top 1500 amplification driven and 1500 deletion driven genes based on fdr

# amplification driven genes
iGCAmp <- read.delim("PancreaticCancer/iGC/iGC_PAAD_GainDriven_genes.tsv",header = T, sep = "\t")
iGCAmp <- iGCAmp[,c(1,3)]
iGCAmp <- iGCAmp[iGCAmp$fdr<0.05,]

# deletion driven genes
iGCDel <- read.delim("PancreaticCancer/iGC/iGC_PAAD_LossDriven_genes.tsv",header = T, sep = "\t")
iGCDel <- iGCDel[,c(1,3)]
iGCDel <- iGCDel[iGCDel$fdr<0.05,]

# Both genes
iGCBoth <- read.delim("PancreaticCancer/iGC/CNA_BothDriven_genes.tsv",header = T, sep = "\t")
iGCBothAmp <- iGCBoth[(iGCBoth$gain_fdr<0.05),]
iGCBothAmp <- iGCBothAmp[,c(1,3)]
colnames(iGCBothAmp)[2] <- "fdr"

iGCBothDel <- iGCBoth[(iGCBoth$loss_fdr<0.05),]
iGCBothDel <- iGCBothDel[,c(1,5)]
colnames(iGCBothDel)[2] <- "fdr"

iGCAmp <- rbind(iGCAmp, iGCBothAmp)
iGCAmpSelected <- iGCAmp # less than 1500
iGCDel <- rbind(iGCDel, iGCBothDel)
iGCDel <- iGCDel[order(iGCDel$fdr),]
iGCDelSelected <- iGCDel[1:1500,]

# Both genes and Total Genes
iGCBoth <- intersect(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE))
iGCTotal <- unique(c(as.character(iGCAmpSelected$GENE), as.character(iGCDelSelected$GENE)))

# save genes
write.table(iGCTotal, file = "PancreaticCancer/iGC/iGC_PAAD_3000Genes.txt", row.names = F, quote = F, sep = "\t")

####################################################################################################################
## plrs
# plrs does not output amplfied and deletion driven genes. Hence selected top 3000 genes from results based on fdr

# total genes
plrsGenes <- read.delim("PancreaticCancer/PLRS/PLRS_PAAD_Results.tsv", header = T, sep = "\t")
plrsGenes<- plrsGenes[plrsGenes$BH.adj.pval<0.05,] 
plrsGenes <- plrsGenes[order(plrsGenes$BH.adj.pval),]
plrsSelected <- plrsGenes[1:3000,]
plrsTotal <- as.character(plrsSelected$Gene)

# save genes
write.table(plrsTotal, file = "PancreaticCancer/PLRS/PLRS_CompPAAD_3000Genes.txt", row.names = F, quote = F, sep = "\t")

####################################################################################################################
## Oncodrive-CIS
# did not provide p-values. Hence we estimate p-values using z-score and estimate fdr
# select top 1500 amplification driven and 1500 deletion driven gene based on fdr

# amplification driven genes
OncoAmp <- read.delim("PancreaticCancer/OncodriveCIS/OncoCNA.AMP",
                      header = T, sep = "\t")
OncoAmp$pval <- pnorm(-abs(OncoAmp$Z.COMBINED))
OncoAmp$p.adj.fdr <- p.adjust(OncoAmp$pval, method = "fdr")
OncoAmp <- OncoAmp[OncoAmp$Z.COMBINED>0,]
OncoAmpSelected <- OncoAmp[OncoAmp$p.adj.fdr<0.05,] # less than 1500.. so using all

# deletion driven genes
OncoDel <- read.delim("PancreaticCancer/OncodriveCIS/OncoCNA.DEL",
                      header = T, sep = "\t")
head(OncoDel)
tail(OncoDel)
OncoDel$pval <- pnorm(-abs(OncoDel$Z.COMBINED))
OncoDel$p.adj.fdr <- p.adjust(OncoDel$pval, method = "fdr")
OncoDel <- OncoDel[OncoDel$Z.COMBINED<0,]
OncoDelSelected <- OncoDel[OncoDel$p.adj.fdr<0.05,] # less than 1500 .. so using all

# Both genes
OncoBoth <- intersect(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID))

# Total genes
OncoTotal <- unique(c(as.character(OncoAmpSelected$ENS.ID), as.character(OncoDelSelected$ENS.ID)))

# save genes
write.table(OncoTotal, file = "PancreaticCancer/OncodriveCIS/Oncodrive_PAAD_3000Genes.txt",
            row.names = F, quote = F, sep = "\t")

####################################################################################################################
## cBioPortal Cancer Gene List and Venn
# copied the list of genes from http://www.cbioportal.org/cancer_gene_list.jsp onto a text file

library(readr)
# reading the text file
mystring <- read_file("PancreaticCancer/Genelevel_Comparison/cBioPortalCancerGeneList.txt")
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
setwd("PancreaticCancer/Genelevel_Comprison")
jpeg("PAAD_3000GenesVenn.jpg", width = 1260, height = 890)  
venn(list(CNAmet = CNAmetTotal,  iGC = iGCTotal,  PLRS = plrsTotal,
          OncodriveCIS = OncoTotal,
          CCGL = cancergenes),
     zcolor = "style", cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()

