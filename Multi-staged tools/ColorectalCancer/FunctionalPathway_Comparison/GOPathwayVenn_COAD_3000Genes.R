# Venns for Gene Ontology BP, CC, MF (DAVID) and Pathway (Reactome) - Colorectal cancer datsaet
# The GO terms for BP, MF, and CC were identified using DAVID 6.8 using the top 1500 amplification 
# and 1500 deletion drive genes
# The pathway terms were identified using Reactome Pathway

library(venn)
rm(list = ls())

## Biological Process (BP) terms Venn analysis
# Reading the files
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/")

# cancer genes BP terms
cancerBP <- read.delim("Venns/DAVID/CancerGenes_GOTERM_BP_FAT.txt",
                       header = T, sep = "\t")
cancerBP <- cancerBP[cancerBP$FDR<0.05,]
cancerBP <- as.character(cancerBP$Term)

# CNAmet genes BP terms
CNAmetBP <- read.delim("CNAmet/DAVID/CNAmet_CompCOAD_3000Genes_GOTERM_BP_FAT.txt",
                       header = T, sep = "\t")
CNAmetBP <- CNAmetBP[CNAmetBP$FDR<0.05,]
CNAmetBP <- as.character(CNAmetBP$Term)

#iGC genes BP terms
iGCBP <- read.delim("iGC/DAVID/iGC_CompCOAD_3000Genes_GOTERM_BP_FAT.txt",
                    header = T, sep = "\t")
iGCBP <- iGCBP[iGCBP$FDR<0.05,]
iGCBP <- as.character(iGCBP$Term)

# PLRS genes BP terms
plrsBP <- read.delim("plrs/DAVID/PLRS_CompCOAD_3000Genes_GOTERM_BP_FAT.txt",
                     header = T, sep = "\t")
plrsBP <- plrsBP[plrsBP$FDR<0.05,]
plrsBP <- as.character(plrsBP$Term)

# Oncodrive-CIS genes BP terms
OncoBP <- read.delim("OncodriveCIS/DAVID/Oncodrive_CompCOAD_3000Genes_GOTERM_BP_FAT.txt",
                     header = T, sep = "\t")
OncoBP <- OncoBP[OncoBP$FDR<0.05,]
OncoBP <- as.character(OncoBP$Term)

## intersections and Venns
cancerCNAmetBP <- intersect(cancerBP, CNAmetBP)
canceriGCBP <- intersect(cancerBP, iGCBP)
cancerplrsBP <- intersect(cancerBP, plrsBP)
cancerOncoBP <- intersect(cancerBP, OncoBP)

# saving the venn diagram
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/Venns/")
jpeg("COAD_3000BPVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetBP,  iGC = iGCBP,  PLRS = plrsBP,
          OncodriveCIS = OncoBP,
          CCGL = cancerBP),
     zcolor = "style", cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()

####################################################################################################
## Cellular Components (CC) terms Venn analysis

setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/")

# Cancer genes
cancerCC <- read.delim("Venns/DAVID/CancerGenes_GOTERM_CC_FAT.txt",
                       header = T, sep = "\t")
cancerCC <- cancerCC[cancerCC$FDR<0.05,]
cancerCC <- as.character(cancerCC$Term)

# CNAmet
CNAmetCC <- read.delim("CNAmet/DAVID/CNAmet_CompCOAD_3000Genes_GOTERM_CC_FAT.txt",
                       header = T, sep = "\t")
CNAmetCC <- CNAmetCC[CNAmetCC$FDR<0.05,]
CNAmetCC <- as.character(CNAmetCC$Term)

# iGC
iGCCC <- read.delim("iGC/DAVID/iGC_CompCOAD_3000Genes_GOTERM_CC_FAT.txt",
                    header = T, sep = "\t")
iGCCC <- iGCCC[iGCCC$FDR<0.05,]
iGCCC <- as.character(iGCCC$Term)

# PLRS
plrsCC <- read.delim("plrs/DAVID/PLRS_CompCOAD_3000Genes_GOTERM_CC_FAT.txt",
                     header = T, sep = "\t")
plrsCC <- plrsCC[plrsCC$FDR<0.05,]
plrsCC <- as.character(plrsCC$Term)

# Oncodrive-CIS
OncoCC <- read.delim("OncodriveCIS/DAVID/Oncodrive_CompCOAD_3000Genes_GOTERM_CC_FAT.txt",
                     header = T, sep = "\t")
OncoCC <- OncoCC[OncoCC$FDR<0.05,]
OncoCC <- as.character(OncoCC$Term)

## intersections and Venns
cancerCNAmetCC <- intersect(cancerCC, CNAmetCC)
canceriGCCC <- intersect(cancerCC, iGCCC)
cancerplrsCC <- intersect(cancerCC, plrsCC)
cancerOncoCC <- intersect(cancerCC, OncoCC)

# saving the Venn diagram
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/Venns/")
jpeg("COAD_3000CCVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetCC,  iGC = iGCCC,  PLRS = plrsCC,
          OncodriveCIS = OncoCC,
          CCGL = cancerCC),
     zcolor = "style",cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()

########################################################################################################
# Molecular Functional (MF) Venn analysis

setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/")

# cancer
cancerMF <- read.delim("Venns/DAVID/CancerGenes_GOTERM_MF_FAT.txt",
                       header = T, sep = "\t")
cancerMF <- cancerMF[cancerMF$FDR<0.05,]
cancerMF <- as.character(cancerMF$Term)

# CNAmet
CNAmetMF <- read.delim("CNAmet/DAVID/CNAmet_CompCOAD_3000Genes_GOTERM_MF_FAT.txt",
                       header = T, sep = "\t")
CNAmetMF <- CNAmetMF[CNAmetMF$FDR<0.05,]
CNAmetMF <- as.character(CNAmetMF$Term)

# iGC
iGCMF <- read.delim("iGC/DAVID/iGC_CompCOAD_3000Genes_GOTERM_MF_FAT.txt",
                    header = T, sep = "\t")
iGCMF <- iGCMF[iGCMF$FDR<0.05,]
iGCMF <- as.character(iGCMF$Term)

# PLRS
plrsMF <- read.delim("plrs/DAVID/PLRS_CompCOAD_3000Genes_GOTERM_MF_FAT.txt",
                     header = T, sep = "\t")
plrsMF <- plrsMF[plrsMF$FDR<0.05,]
plrsMF <- as.character(plrsMF$Term)

# Oncodrive-CIS
OncoMF <- read.delim("OncodriveCIS/DAVID/Oncodrive_CompCOAD_3000Genes_GOTERM_MF_FAT.txt",
                     header = T, sep = "\t")
OncoMF <- OncoMF[OncoMF$FDR<0.05,]
OncoMF <- as.character(OncoMF$Term)

## intersections and Venns
cancerCNAmetMF <- intersect(cancerMF, CNAmetMF)
canceriGCMF <- intersect(cancerMF, iGCMF)
cancerplrsMF <- intersect(cancerMF, plrsMF)
cancerOncoMF <- intersect(cancerMF, OncoMF)

# saving venn diagrams
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/Venns/")
jpeg("COAD_3000MFVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetMF,  iGC = iGCMF,  PLRS = plrsMF,
          OncodriveCIS = OncoMF,
          CCGL = cancerMF),
     zcolor = "style", cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()

####################################################################################################
# Pathways Venn Analysis
# Reading the files
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/")
# Cancer
cancerPathway <- read.delim("Venns/DAVID/CancerGenes_Reactome_Pathway.csv",
                            header = T, sep = ",")
cancerPathway <- cancerPathway[cancerPathway$Entities.FDR<0.05,]
cancerPathway <- as.character(cancerPathway$Pathway.identifier)

#CNAmet
CNAmetPathway <- read.delim("CNAmet/DAVID/CNAmet_CompCOAD_3000Genes_Reactome_Pathway.csv",
                            header = T, sep = ",")
CNAmetPathway <- CNAmetPathway[CNAmetPathway$Entities.FDR<0.05,]
CNAmetPathway <- as.character(CNAmetPathway$Pathway.identifier)

# iGC
iGCPathway <- read.delim("iGC/DAVID/iGC_CompCOAD_3000Genes_Reactome_Pathway.csv",
                         header = T, sep = ",")
iGCPathway <- iGCPathway[iGCPathway$Entities.FDR<0.05,]
iGCPathway <- as.character(iGCPathway$Pathway.identifier)

# PLRS
plrsPathway <- read.delim("plrs/DAVID/PLRS_CompCOAD_3000Genes_Reactome_Pathway.csv",
                          header = T, sep = ",")
plrsPathway <- plrsPathway[plrsPathway$Entities.FDR<0.05,]
plrsPathway <- as.character(plrsPathway$Pathway.identifier)

# Oncodrive-CIS
OncoPathway <- read.delim("OncodriveCIS/DAVID/Oncodrive_CompCOAD_3000Genes_Reactome_Pathway.csv",
                          header = T, sep = ",")
OncoPathway <- OncoPathway[OncoPathway$Entities.FDR<0.05,]
OncoPathway <- as.character(OncoPathway$Pathway.identifier)

## intersections and Venns
cancerCNAmetPathway <- intersect(cancerPathway, CNAmetPathway)
canceriGCPathway <- intersect(cancerPathway, iGCPathway)
cancerplrsPathway <- intersect(cancerPathway, plrsPathway)
cancerOncoPathway <- intersect(cancerPathway, OncoPathway)

# Saving the venn diagram
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/Venns/")
jpeg("COAD_3000PathwayVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetPathway,  iGC = iGCPathway,  PLRS = plrsPathway,
          OncodriveCIS = OncoPathway,
          CCGL = cancerPathway),
     zcolor = "style", cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()
