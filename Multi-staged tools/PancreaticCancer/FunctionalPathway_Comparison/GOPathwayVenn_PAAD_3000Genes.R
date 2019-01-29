# Venns for Gene Ontology BP, CC, MF (DAVID) and Pathway (Reactome) - Colorectal cancer datsaet
# The GO terms for BP, MF, and CC were identified using DAVID 6.8 using the top 1500 amplification 
# and 1500 deletion drive genes
# The pathway terms were identified using Reactome Pathway

library(venn)
rm(list = ls())

## Biological Process (BP) terms Venn analysis
# Reading the files

# cancer genes BP terms
cancerBP <- read.delim("FunctionalPathway_Comparison/DAVID/CancerGenes_GOTERM_BP_FAT.txt",
                       header = T, sep = "\t")
cancerBP <- cancerBP[cancerBP$FDR<0.05,]
cancerBP <- as.character(cancerBP$Term)

# CNAmet genes BP terms
CNAmetBP <- read.delim("PancreaticCancer/CNAmet/DAVID/CNAmet_PAAD_3000Genes_GOTERM_BP_FAT.txt",
                       header = T, sep = "\t")
CNAmetBP <- CNAmetBP[CNAmetBP$FDR<0.05,]
CNAmetBP <- as.character(CNAmetBP$Term)

#iGC genes BP terms
iGCBP <- read.delim("PancreaticCancer/iGC/DAVID/iGC_PAAD_3000Genes_GOTERM_BP_FAT.txt",
                    header = T, sep = "\t")
iGCBP <- iGCBP[iGCBP$FDR<0.05,]
iGCBP <- as.character(iGCBP$Term)

# PLRS genes BP terms
plrsBP <- read.delim("PancreaticCancer/PLRS/DAVID/PLRS_PAAD_3000Genes_GOTERM_BP_FAT.txt",
                     header = T, sep = "\t")
plrsBP <- plrsBP[plrsBP$FDR<0.05,]
plrsBP <- as.character(plrsBP$Term)

# Oncodrive-CIS genes BP terms
OncoBP <- read.delim("PancreaticCancer/OncodriveCIS/DAVID/Oncodrive_PAAD_3000Genes_GOTERM_BP_FAT.txt",
                     header = T, sep = "\t")
OncoBP <- OncoBP[OncoBP$FDR<0.05,]
OncoBP <- as.character(OncoBP$Term)

## intersections and Venns
cancerCNAmetBP <- intersect(cancerBP, CNAmetBP)
canceriGCBP <- intersect(cancerBP, iGCBP)
cancerplrsBP <- intersect(cancerBP, plrsBP)
cancerOncoBP <- intersect(cancerBP, OncoBP)

# saving the venn diagram
jpeg("PAAD_3000BPVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetBP,  iGC = iGCBP,  PLRS = plrsBP,
          OncodriveCIS = OncoBP,
          CCGL = cancerBP),
     zcolor = "style", cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()

####################################################################################################
## Cellular Components (CC) terms Venn analysis

# cancer genes
cancerCC <- read.delim("FunctionalPathway_Comparison/DAVID/CancerGenes_GOTERM_CC_FAT.txt",
                       header = T, sep = "\t")
cancerCC <- cancerCC[cancerCC$FDR<0.05,]
cancerCC <- as.character(cancerCC$Term)

# CNAmet
CNAmetCC <- read.delim("PancreaticCancer/CNAmet/DAVID/CNAmet_PAAD_3000Genes_GOTERM_CC_FAT.txt",
                       header = T, sep = "\t")
CNAmetCC <- CNAmetCC[CNAmetCC$FDR<0.05,]
CNAmetCC <- as.character(CNAmetCC$Term)

#iGC
iGCCC <- read.delim("PancreaticCancer/iGC/DAVID/iGC_PAAD_3000Genes_GOTERM_CC_FAT.txt",
                    header = T, sep = "\t")
iGCCC <- iGCCC[iGCCC$FDR<0.05,]
iGCCC <- as.character(iGCCC$Term)

# PLRS
plrsCC <- read.delim("PancreaticCancer/PLRS/DAVID/PLRS_PAAD_3000Genes_GOTERM_CC_FAT.txt",
                     header = T, sep = "\t")
plrsCC <- plrsCC[plrsCC$FDR<0.05,]
plrsCC <- as.character(plrsCC$Term)

# Oncodrive-CIS
OncoCC <- read.delim("PancreaticCancer/OncodriveCIS/DAVID/Oncodrive_PAAD_3000Genes_GOTERM_CC_FAT.txt",
                     header = T, sep = "\t")
OncoCC <- OncoCC[OncoCC$FDR<0.05,]
OncoCC <- as.character(OncoCC$Term)

## intersections and Venns
cancerCNAmetCC <- intersect(cancerCC, CNAmetCC)
canceriGCCC <- intersect(cancerCC, iGCCC)
cancerplrsCC <- intersect(cancerCC, plrsCC)
cancerOncoCC <- intersect(cancerCC, OncoCC)

# Saving the Venn diagram
jpeg("PAAD_3000CCVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetCC,  iGC = iGCCC,  PLRS = plrsCC,
          OncodriveCIS = OncoCC,
          CCGL = cancerCC),
     zcolor = "style",cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()

########################################################################################################
# Molecular Functional (MF) Venn analysis

# cancer
cancerMF <- read.delim("FunctionalPathway_Comparison/DAVID/CancerGenes_GOTERM_MF_FAT.txt",
                       header = T, sep = "\t")
cancerMF <- cancerMF[cancerMF$FDR<0.05,]
cancerMF <- as.character(cancerMF$Term)

# CNAmet
CNAmetMF <- read.delim("PancreaticCancer/CNAmet/DAVID/CNAmet_PAAD_3000Genes_GOTERM_MF_FAT.txt",
                       header = T, sep = "\t")
CNAmetMF <- CNAmetMF[CNAmetMF$FDR<0.05,]
CNAmetMF <- as.character(CNAmetMF$Term)

# iGC
iGCMF <- read.delim("PancreaticCancer/iGC/DAVID/iGC_PAAD_3000Genes_GOTERM_MF_FAT.txt",
                    header = T, sep = "\t")
iGCMF <- iGCMF[iGCMF$FDR<0.05,]
iGCMF <- as.character(iGCMF$Term)

# PLRS
plrsMF <- read.delim("PancreaticCancer/PLRS/DAVID/PLRS_PAAD_3000Genes_GOTERM_MF_FAT.txt",
                     header = T, sep = "\t")
plrsMF <- plrsMF[plrsMF$FDR<0.05,]
plrsMF <- as.character(plrsMF$Term)

# Oncodrive-CIS
OncoMF <- read.delim("PancreaticCancer/OncodriveCIS/DAVID/Oncodrive_CompPAAD_3000Genes_GOTERM_MF_FAT.txt",
                     header = T, sep = "\t")
OncoMF <- OncoMF[OncoMF$FDR<0.05,]
OncoMF <- as.character(OncoMF$Term)

## intersections and Venns
cancerCNAmetMF <- intersect(cancerMF, CNAmetMF)
canceriGCMF <- intersect(cancerMF, iGCMF)
cancerplrsMF <- intersect(cancerMF, plrsMF)
cancerOncoMF <- intersect(cancerMF, OncoMF)


# saving the Venn
jpeg("PAAD_3000MFVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetMF,  iGC = iGCMF,  PLRS = plrsMF,
          OncodriveCIS = OncoMF,
          CCGL = cancerMF),
     zcolor = "style", cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()

####################################################################################################
# Pathways Venn Analysis
# Reading the files

# cancer
cancerPathway <- read.delim("FunctionalPathway_Comparison/DAVID/CancerGenes_Reactome_Pathway.csv",
                            header = T, sep = ",")
cancerPathway <- cancerPathway[cancerPathway$Entities.FDR<0.05,]
cancerPathway <- as.character(cancerPathway$Pathway.identifier)

# CNAmet
CNAmetPathway <- read.delim("PancreaticCancer/CNAmet/DAVID/CNAmet_PAAD_3000Genes_Reactome_Pathway.csv",
                            header = T, sep = ",")
CNAmetPathway <- CNAmetPathway[CNAmetPathway$Entities.FDR<0.05,]
CNAmetPathway <- as.character(CNAmetPathway$Pathway.identifier)

# iGC
iGCPathway <- read.delim("PancreaticCancer/iGC/DAVID/iGC_PAAD_3000Genes_Reactome_Pathway.csv",
                         header = T, sep = ",")
iGCPathway <- iGCPathway[iGCPathway$Entities.FDR<0.05,]
iGCPathway <- as.character(iGCPathway$Pathway.identifier)

# PLRS
plrsPathway <- read.delim("PancreaticCancer/PLRS/DAVID/PLRS_PAAD_3000Genes_Reactome_Pathway.csv",
                          header = T, sep = ",")
plrsPathway <- plrsPathway[plrsPathway$Entities.FDR<0.05,]
plrsPathway <- as.character(plrsPathway$Pathway.identifier)

# Oncodrive-CIS
OncoPathway <- read.delim("PancreaticCancer/OncodriveCIS/DAVID/Oncodrive_PAAD_3000Genes_Reactome_Pathway.csv",
                          header = T, sep = ",")
OncoPathway <- OncoPathway[OncoPathway$Entities.FDR<0.05,]
OncoPathway <- as.character(OncoPathway$Pathway.identifier)

## intersections and Venns
cancerCNAmetPathway <- intersect(cancerPathway, CNAmetPathway)
canceriGCPathway <- intersect(cancerPathway, iGCPathway)
cancerplrsPathway <- intersect(cancerPathway, plrsPathway)
cancerOncoPathway <- intersect(cancerPathway, OncoPathway)

# Saving the venn diagrame=
jpeg("PAAD_3000PathwayVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetPathway,  iGC = iGCPathway,  PLRS = plrsPathway,
          OncodriveCIS = OncoPathway,
          CCGL = cancerPathway),
     zcolor = "style", cexil = 2, cexsn = 2.5, opacity = 0.3,size = 20)
dev.off()
