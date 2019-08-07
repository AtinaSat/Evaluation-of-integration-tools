# Venns for Genes, Gene Ontology BP, CC, MF (using clusterProfiler) and Pathway (using Reactome DB) - Pancreatic cancer datsaet
# The GO terms were identified using enrichGO from clusterProfiler
# The reactome pathways were identified using submitting the significant 3000 genes from each tool in Reactome pathway database
# date: 25/06/2019

rm(list = ls())  

library(venn)
library(clusterProfiler)
library(ReactomePA)
source("/home/anita/Benchmarking/two_omics/HUGOtoEntrezIDsFunction.R")

## cBioPortal cancer genes list (CCGL)
# Reading the file
library(readr)
mystring <- read_file("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns/cBioPortalCancerGeneList.txt")
cBioTotal <- strsplit(mystring, "\n", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
cBioTotal[cBioTotal=="CDKN2AP14ARF"] <- "CDKN2AP14ARF"
cBioEntrez <- HUGOtoEntrez(InputGenes = cBioTotal)
# GO enrichment
cBioGO <- enrichGO(gene = cBioEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
cBioGOres <- cBioGO@result
cBioBP <- cBioGOres[cBioGOres$ONTOLOGY=="BP",]
cBioCC <- cBioGOres[cBioGOres$ONTOLOGY=="CC",]
cBioMF <- cBioGOres[cBioGOres$ONTOLOGY=="MF",]
# Pathway enrichment
cBioPathway <- enrichPathway(gene = cBioEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

EstimateCCGLOverlap <- function(Genes, BP, CC, MF, Pathway){
  CommonGenes <- length(intersect(Genes, cBioTotal))
  CommonBP <- length(intersect(BP$ID, cBioBP$ID))
  CommonCC <- length(intersect(CC$ID, cBioCC$ID))
  CommonMF <- length(intersect(MF$ID, cBioMF$ID))
  CommonPathway <- length(intersect(Pathway, cBioPathway))
  return(c(CommonGenes, CommonBP, CommonCC, CommonMF, CommonPathway))
}


## Cancer Gene Consensus genes
# Reading the file
CGCTotal <- read.csv("/home/anita/Benchmarking/two_omics/Census_allThu Jun 20 07_06_05 2019.csv", header = T)
CGCTotal <- as.character(CGCTotal$Gene.Symbol)
CGCEntrez <- HUGOtoEntrez(InputGenes = CGCTotal)
# GO enrichment
CGCGO <- enrichGO(gene = CGCEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
CGCGOres <- CGCGO@result
CGCBP <- CGCGOres[CGCGOres$ONTOLOGY=="BP",]
CGCCC <- CGCGOres[CGCGOres$ONTOLOGY=="CC",]
CGCMF <- CGCGOres[CGCGOres$ONTOLOGY=="MF",]
# Pathway enrichment
CGCPathway <- enrichPathway(gene = CGCEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

EstimateCGCOverlap <- function(Genes, BP, CC, MF, Pathway){
  CommonGenes <- length(intersect(Genes, CGCTotal))
  CommonBP <- length(intersect(BP$ID, CGCBP$ID))
  CommonCC <- length(intersect(CC$ID, CGCCC$ID))
  CommonMF <- length(intersect(MF$ID, CGCMF$ID))
  CommonPathway <- length(intersect(Pathway, CGCPathway))
  return(c(CommonGenes, CommonBP, CommonCC, CommonMF, CommonPathway))
}

## creating a dataframe to fill in the results as we progress through the script
PAAD_FunctionalAnalysis <- data.frame(matrix(ncol = 17, nrow = 6))
colnames(PAAD_FunctionalAnalysis) <- c("Tools", "Omics", "Sig.Genes", "Sig.Genes in CCGL", "Sig.Genes in CGC",
                                       "Sig.BP", "Sig.BP in CCGL BP", "Sig.BP in CGC BP",
                                       "Sig.CC", "Sig.BP in CCGL CC", "Sig.BP in CGC CC",
                                       "Sig.MF", "Sig.BP in CCGL MF", "Sig.BP in CGC MF",
                                       "Sig.Pathways", "Sig.Pathways in CCGL Pathways", "Sig.Pathways in CGC Pathways")

PAAD_FunctionalAnalysis$Tools <- c("CNAMet", "iGC", "PLRS", "Oncodrive-CIS", "CNAMet", "MethylMix")
PAAD_FunctionalAnalysis$Omics <- c(rep("CNA+GE", 4), rep("ME+GE", 2))

## CNAmet
# Reading the file
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/")
CNAmetTotal <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/CNAmet_PAAD_3000Genes_OncodriveInput.txt",
                          header = T, sep = "\t") 
CNAmetTotal <- as.character(CNAmetTotal$x)
CNAmetEntrez <- HUGOtoEntrez(InputGenes = CNAmetTotal)
# GO enrichment
CNAmetGO <- enrichGO(gene = CNAmetEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
CNAmetGOres <- CNAmetGO@result
CNAmetBP <- CNAmetGOres[CNAmetGOres$ONTOLOGY=="BP",]
CNAmetCC <- CNAmetGOres[CNAmetGOres$ONTOLOGY=="CC",]
CNAmetMF <- CNAmetGOres[CNAmetGOres$ONTOLOGY=="MF",]
# Pathway enrichment
CNAmetPathway <- enrichPathway(gene = CNAmetEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

PAAD_FunctionalAnalysis[1, c(3, 6, 9, 12, 15)] <- c(length(CNAmetTotal), nrow(CNAmetBP), nrow(CNAmetCC), nrow(CNAmetMF), length(CNAmetPathway)) # Significant genes
PAAD_FunctionalAnalysis[1, c(4, 7, 10, 13, 16)] <- EstimateCCGLOverlap(Genes = CNAmetTotal, BP = CNAmetBP, CC = CNAmetCC, MF = CNAmetMF,
                                                                       Pathway = CNAmetPathway) # Overlap with CCGL
PAAD_FunctionalAnalysis[1, c(5, 8, 11, 14, 17)] <- EstimateCGCOverlap(Genes = CNAmetTotal, BP = CNAmetBP, CC = CNAmetCC, MF = CNAmetMF,
                                                                      Pathway = CNAmetPathway) # Overlap with CGC

## iGC
# Reading the file
iGCTotal <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/iGC_PAAD_3000Genes_OncodriveInput.txt", header = T, sep = "\t")
iGCTotal <- as.character(iGCTotal$x)
iGCEntrez <- HUGOtoEntrez(InputGenes = iGCTotal)
# GO enrichment
iGCGO <- enrichGO(gene = iGCEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
iGCGOres <- iGCGO@result
iGCBP <- iGCGOres[iGCGOres$ONTOLOGY=="BP",]
iGCCC <- iGCGOres[iGCGOres$ONTOLOGY=="CC",]
iGCMF <- iGCGOres[iGCGOres$ONTOLOGY=="MF",]
# Pathway enrichment
iGCPathway <- enrichPathway(gene = iGCEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

PAAD_FunctionalAnalysis[2, c(3, 6, 9, 12, 15)] <- c(length(iGCTotal), nrow(iGCBP), nrow(iGCCC), nrow(iGCMF), length(iGCPathway)) # Significant genes
PAAD_FunctionalAnalysis[2, c(4, 7, 10, 13, 16)] <- EstimateCCGLOverlap(Genes = iGCTotal, BP = iGCBP, CC = iGCCC, MF = iGCMF,
                                                                       Pathway = iGCPathway) # Overlap with CCGL
PAAD_FunctionalAnalysis[2, c(5, 8, 11, 14, 17)] <- EstimateCGCOverlap(Genes = iGCTotal, BP = iGCBP, CC = iGCCC, MF = iGCMF,
                                                                      Pathway = iGCPathway) # Overlap with CGC


## PLRS 
# Reading the file
plrsTotal <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/PLRS_PAAD_3000Genes_OncodriveInput.txt", header = T, sep = "\t")
plrsTotal <- as.character(plrsTotal$x)
plrsEntrez <- HUGOtoEntrez(InputGenes = plrsTotal)
# GO enrichment
plrsGO <- enrichGO(gene = plrsEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
plrsGOres <- plrsGO@result
plrsBP <- plrsGOres[plrsGOres$ONTOLOGY=="BP",]
plrsCC <- plrsGOres[plrsGOres$ONTOLOGY=="CC",]
plrsMF <- plrsGOres[plrsGOres$ONTOLOGY=="MF",]
# Pathway enrichment
plrsPathway <- enrichPathway(gene = plrsEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

PAAD_FunctionalAnalysis[3, c(3, 6, 9, 12, 15)] <- c(length(plrsTotal), nrow(plrsBP), nrow(plrsCC), nrow(plrsMF), length(plrsPathway)) # Significant genes
PAAD_FunctionalAnalysis[3, c(4, 7, 10, 13, 16)] <- EstimateCCGLOverlap(Genes = plrsTotal, BP = plrsBP, CC = plrsCC, MF = plrsMF,
                                                                       Pathway = plrsPathway) # Overlap with CCGL
PAAD_FunctionalAnalysis[3, c(5, 8, 11, 14, 17)] <- EstimateCGCOverlap(Genes = plrsTotal, BP = plrsBP, CC = plrsCC, MF = plrsMF,
                                                                      Pathway = plrsPathway) # Overlap with CGC

## Oncodrive-CIS 
# Reading the file
OncoTotal <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/OncodriveCIS/Oncodrive_PAAD_3000Genes.txt", header = T, sep = "\t")
OncoTotal <- as.character(OncoTotal$x)
OncoEntrez <- HUGOtoEntrez(InputGenes = OncoTotal)
# GO enrichment
OncoGO <- enrichGO(gene = OncoEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
OncoGOres <- OncoGO@result
OncoBP <- OncoGOres[OncoGOres$ONTOLOGY=="BP",]
OncoCC <- OncoGOres[OncoGOres$ONTOLOGY=="CC",]
OncoMF <- OncoGOres[OncoGOres$ONTOLOGY=="MF",]
# Pathway enrichment
OncoPathway <- enrichPathway(gene = OncoEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

PAAD_FunctionalAnalysis[4, c(3, 6, 9, 12, 15)] <- c(length(OncoTotal), nrow(OncoBP), nrow(OncoCC), nrow(OncoMF), length(OncoPathway)) # Significant genes
PAAD_FunctionalAnalysis[4, c(4, 7, 10, 13, 16)] <- EstimateCCGLOverlap(Genes = OncoTotal, BP = OncoBP, CC = OncoCC, MF = OncoMF,
                                                                       Pathway = OncoPathway) # Overlap with CCGL
PAAD_FunctionalAnalysis[4, c(5, 8, 11, 14, 17)] <- EstimateCGCOverlap(Genes = OncoTotal, BP = OncoBP, CC = OncoCC, MF = OncoMF,
                                                                      Pathway = OncoPathway) # Overlap with CGC


#################################################################################################################################
##### Methylation Results
## CNAmet - Methyl
# Reading the file
CNAmetMethylTotal <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/CNAmetMethyl/CNAmetMethyl_PAAD_3000Genes.txt",
                                header = T, sep = "\t") 
CNAmetMethylTotal <- as.character(CNAmetMethylTotal$x)
CNAmetMethylEntrez <- HUGOtoEntrez(InputGenes = CNAmetMethylTotal)
# GO enrichment
CNAmetMethylGO <- enrichGO(gene = CNAmetMethylEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
CNAmetMethylGOres <- CNAmetMethylGO@result
CNAmetMethylBP <- CNAmetMethylGOres[CNAmetMethylGOres$ONTOLOGY=="BP",]
CNAmetMethylCC <- CNAmetMethylGOres[CNAmetMethylGOres$ONTOLOGY=="CC",]
CNAmetMethylMF <- CNAmetMethylGOres[CNAmetMethylGOres$ONTOLOGY=="MF",]
# Pathway enrichmet
CNAmetMethylPathway <- enrichPathway(gene = CNAmetMethylEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

PAAD_FunctionalAnalysis[5, c(3, 6, 9, 12, 15)] <- c(length(CNAmetMethylTotal), nrow(CNAmetMethylBP), nrow(CNAmetMethylCC), nrow(CNAmetMethylMF), length(CNAmetMethylPathway)) # Significant genes
PAAD_FunctionalAnalysis[5, c(4, 7, 10, 13, 16)] <- EstimateCCGLOverlap(Genes = CNAmetMethylTotal, BP = CNAmetMethylBP, CC = CNAmetMethylCC, MF = CNAmetMethylMF,
                                                                       Pathway = CNAmetMethylPathway) # Overlap with CCGL
PAAD_FunctionalAnalysis[5, c(5, 8, 11, 14, 17)] <- EstimateCGCOverlap(Genes = CNAmetMethylTotal, BP = CNAmetMethylBP, CC = CNAmetMethylCC, MF = CNAmetMethylMF,
                                                                      Pathway = CNAmetMethylPathway) # Overlap with CGC


## MethylMix
# Reading the file
MethylMixTotal <- read.delim("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/MethylMix/MethylMix_PAAD_3000Genes.txt",
                             header = T, sep = "\t") 
MethylMixTotal <- as.character(MethylMixTotal$x)
MethylMixEntrez <- HUGOtoEntrez(InputGenes = MethylMixTotal)
# GO enrichment
MethylMixGO <- enrichGO(gene = MethylMixEntrez$Entrez, OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
MethylMixGOres <- MethylMixGO@result
MethylMixBP <- MethylMixGOres[MethylMixGOres$ONTOLOGY=="BP",]
MethylMixCC <- MethylMixGOres[MethylMixGOres$ONTOLOGY=="CC",]
MethylMixMF <- MethylMixGOres[MethylMixGOres$ONTOLOGY=="MF",]
# Pathway enrichment
MethylMixPathway <- enrichPathway(gene = MethylMixEntrez$Entrez, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr")@result$ID

PAAD_FunctionalAnalysis[6, c(3, 6, 9, 12, 15)] <- c(length(MethylMixTotal), nrow(MethylMixBP), nrow(MethylMixCC), nrow(MethylMixMF), length(MethylMixPathway)) # Significant genes
PAAD_FunctionalAnalysis[6, c(4, 7, 10, 13, 16)] <- EstimateCCGLOverlap(Genes = MethylMixTotal, BP = MethylMixBP, CC = MethylMixCC, MF = MethylMixMF,
                                                                       Pathway = MethylMixPathway) # Overlap with CCGL
PAAD_FunctionalAnalysis[6, c(5, 8, 11, 14, 17)] <- EstimateCGCOverlap(Genes = MethylMixTotal, BP = MethylMixBP, CC = MethylMixCC, MF = MethylMixMF,
                                                                      Pathway = MethylMixPathway) # Overlap with CGC


setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns_OncodriveInput/")
write.table(x = PAAD_FunctionalAnalysis, file = "PAAD_FunctionalAnalysis_withCGC_OncodriveInput_3000Genes_RPackages.txt", sep = "\t", quote = F, row.names = FALSE)
PAAD_FunctionalAnalysis2 <- PAAD_FunctionalAnalysis[, c(1,2,3,4,6,7,9,10,12,13,15,16)] # extracting information relavent only to CCGL
write.table(x = PAAD_FunctionalAnalysis2, file = "PAAD_FunctionalAnalysis_OncodriveInput_3000Genes_RPackages.txt", sep = "\t", quote = F, row.names = FALSE)
save.image(file = "PAAD_FunctionalAnalysis_Results_OncodriveInput_3000Genes_RPackages.RData")

## Venn Plots
# Gene Plot
tiff("PAAD_OncodriveInput_GenesVenn.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetTotal, iGC = iGCTotal, PLRS = plrsTotal, OncodriveCIS = OncoTotal, CCGL = cBioTotal),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# BP plot
tiff("PAAD_OncodriveInput_BPVenn.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetBP$ID, iGC = iGCBP$ID, PLRS = plrsBP$ID, OncodriveCIS = OncoBP$ID, CCGL = cBioBP$ID),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# CC plot
tiff("PAAD_OncodriveInput_CCVenn.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetCC$ID, iGC = iGCCC$ID, PLRS = plrsCC$ID, OncodriveCIS = OncoCC$ID, CCGL = cBioCC$ID),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# MF plot
tiff("PAAD_OncodriveInput_MFVenn.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetMF$ID, iGC = iGCMF$ID, PLRS = plrsMF$ID, OncodriveCIS = OncoMF$ID, CCGL = cBioMF$ID),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# Pathway plot
tiff("PAAD_OncodriveInput_pathwayVenn.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetPathway, iGC = iGCPathway, PLRS = plrsPathway, OncodriveCIS = OncoPathway, CCGL = cBioPathway),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()

length(intersect(intersect(intersect(plrsTotal, iGCTotal), OncoTotal), CNAmetTotal))

length(intersect(plrsTotal, iGCTotal)) # 1630
length(intersect(plrsTotal, OncoTotal)) # 1227
length(intersect(plrsTotal, CNAmetTotal)) # 1635

length(intersect(iGCTotal, OncoTotal)) # 1227
length(intersect(CNAmetTotal, OncoTotal)) # 1270

length(intersect(intersect(CNAmetTotal, plrsTotal), iGCTotal)) # 986
length(intersect(intersect(CNAmetTotal, plrsTotal), OncoTotal)) # 978
length(intersect(intersect(OncoTotal, plrsTotal), iGCTotal)) # 1056
round((length(intersect(intersect(OncoTotal, plrsTotal), iGCTotal))/c(length(iGCTotal), length(OncoTotal), length(plrsTotal)))*100,2)
# [1] 42.32 49.12 35.20



############# Methylation Plots
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns_OncodriveInput/")
load(file = "PAAD_FunctionalAnalysis_Results_OncodriveInput_3000Genes_RPackages.RData")



# Gene Plot
tiff("PAAD_OncodriveInput_GenesVenn_Meth.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetMethylTotal, MethylMix = MethylMixTotal, CCGL = cBioTotal),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# BP plot
tiff("PAAD_OncodriveInput_BPVenn_Meth.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetMethylBP$ID, MethylMix = MethylMixBP$ID, CCGL = cBioBP$ID),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# CC plot
tiff("PAAD_OncodriveInput_CCVenn_Meth.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetMethylCC$ID, MethylMix = MethylMixCC$ID, CCGL = cBioCC$ID),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# MF plot
tiff("PAAD_OncodriveInput_MFVenn_Meth.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetMethylMF$ID, MethylMix = MethylMixMF$ID, CCGL = cBioMF$ID),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()
# Pathway plot
tiff("PAAD_OncodriveInput_pathwayVenn_Meth.tiff", width = 952, height = 952, res = 300)  
venn(x = list(CNAmet = CNAmetMethylPathway, MethylMix = MethylMixPathway, CCGL = cBioPathway),
     zcolor = "style", cexil = 0.6, cexsn = 0.6, opacity = 0.3)
dev.off()




length(cBioBP$ID)
length(CGCBP$ID)
length(CNAmetBP$ID)
length(iGCBP$ID)
length(plrsBP$ID)
length(OncoBP$ID)

length(cBioCC$ID)
length(CGCCC$ID)
length(CNAmetCC$ID)
length(iGCCC$ID)
length(plrsCC$ID)
length(OncoCC$ID)

length(cBioMF$ID)
length(CGCMF$ID)
length(CNAmetMF$ID)
length(iGCMF$ID)
length(plrsMF$ID)
length(OncoMF$ID)


length(intersect(CNAmetBP$ID, cBioBP$ID))
length(intersect(iGCBP$ID, cBioBP$ID))
length(intersect(plrsBP$ID, cBioBP$ID))
length(intersect(OncoBP$ID, cBioBP$ID))

length(intersect(CNAmetCC$ID, cBioCC$ID))
length(intersect(iGCCC$ID, cBioCC$ID))
length(intersect(plrsCC$ID, cBioCC$ID))
length(intersect(OncoCC$ID, cBioCC$ID))

length(intersect(CNAmetMF$ID, cBioMF$ID))
length(intersect(iGCMF$ID, cBioMF$ID))
length(intersect(plrsMF$ID, cBioMF$ID))
length(intersect(OncoMF$ID, cBioMF$ID))


length(intersect(CNAmetBP$ID, CGCBP$ID))
length(intersect(iGCBP$ID, CGCBP$ID))
length(intersect(plrsBP$ID, CGCBP$ID))
length(intersect(OncoBP$ID, CGCBP$ID))

length(intersect(CNAmetCC$ID, CGCCC$ID))
length(intersect(iGCCC$ID, CGCCC$ID))
length(intersect(plrsCC$ID, CGCCC$ID))
length(intersect(OncoCC$ID, CGCCC$ID))

length(intersect(CNAmetMF$ID, CGCMF$ID))
length(intersect(iGCMF$ID, CGCMF$ID))
length(intersect(plrsMF$ID, CGCMF$ID))
length(intersect(OncoMF$ID, CGCMF$ID))

