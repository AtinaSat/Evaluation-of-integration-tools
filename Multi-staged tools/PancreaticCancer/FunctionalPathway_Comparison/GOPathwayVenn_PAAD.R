## Venns for Gene Ontology BP, CC, MF and Pathway
library(venn)
rm(list = ls())

## Biological BP analysis
# Reading the files
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/")
cancerBP <- read.delim("Venns/CancerGenes_GOBP/enrichment_results_wg_result1545723067.txt",
                       header = T, sep = "\t")
cancerBP <- as.character(cancerBP$geneset)

CNAmetBP <- read.delim("CNAmet/CNAmet_compPAAD_GOBP/enrichment_results_wg_result1545719238.txt",
                       header = T, sep = "\t")
CNAmetBP <- as.character(CNAmetBP$geneset)

iGCBP <- read.delim("iGC/iGC_compPAAD_GOBP/enrichment_results_wg_result1545721530.txt",
                    header = T, sep = "\t")
iGCBP <- as.character(iGCBP$geneset)

plrsBP <- read.delim("plrs/PLRS_compPAAD_GOBP/enrichment_results_wg_result1545721979.txt",
                     header = T, sep = "\t")
plrsBP <- as.character(plrsBP$geneset)

OncoBP <- read.delim("OncodriveCIS/Oncodrive_compPAAD_GOBP/enrichment_results_wg_result1545721772.txt",
                     header = T, sep = "\t")
OncoBP <- as.character(OncoBP$geneset)

## intersections and Venns
cancerCNAmetBP <- intersect(cancerBP, CNAmetBP)
canceriGCBP <- intersect(cancerBP, iGCBP)
cancerplrsBP <- intersect(cancerBP, plrsBP)
cancerOncoBP <- intersect(cancerBP, OncoBP)

setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns/")
jpeg("PAAD_BPVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetBP,  iGC = iGCBP,  PLRS = plrsBP,
          OncodriveCIS = OncoBP,
          CCGL = cancerBP),
     zcolor = "style", cexil = 1.5, cexsn = 1.5, opacity = 0.3,size = 20)
dev.off()

####################################################################################################
## Cellular Components Analysis Venn

setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/")
cancerCC <- read.delim("Venns/CancerGenes_GOCC/enrichment_results_wg_result1545723075.txt",
                       header = T, sep = "\t")
cancerCC <- as.character(cancerCC$geneset)

CNAmetCC <- read.delim("CNAmet/CNAmet_compPAAD_GOCC/enrichment_results_wg_result1545719420.txt",
                       header = T, sep = "\t")
CNAmetCC <- as.character(CNAmetCC$geneset)

iGCCC <- read.delim("iGC/iGC_compPAAD_GOCC/enrichment_results_wg_result1545721611.txt",
                    header = T, sep = "\t")
iGCCC <- as.character(iGCCC$geneset)

plrsCC <- read.delim("plrs/PLRS_compPAAD_GOCC/enrichment_results_wg_result1545721990.txt",
                     header = T, sep = "\t")
plrsCC <- as.character(plrsCC$geneset)

OncoCC <- read.delim("OncodriveCIS/Oncodrive_compPAAD_GOCC/enrichment_results_wg_result1545721786.txt",
                     header = T, sep = "\t")
OncoCC <- as.character(OncoCC$geneset)

## intersections and Venns
cancerCNAmetCC <- intersect(cancerCC, CNAmetCC)
canceriGCCC <- intersect(cancerCC, iGCCC)
cancerplrsCC <- intersect(cancerCC, plrsCC)
cancerOncoCC <- intersect(cancerCC, OncoCC)

setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns/")
jpeg("PAAD_CCVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetCC,  iGC = iGCCC,  PLRS = plrsCC,
          OncodriveCIS = OncoCC,
          CCGL = cancerCC),
     zcolor = "style", cexil = 1.5, cexsn = 1.5, opacity = 0.3,size = 20)
dev.off()

########################################################################################################
# Molecular Functional analysis
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/")
cancerMF <- read.delim("Venns/CancerGenes_GOMF/enrichment_results_wg_result1545723089.txt",
                       header = T, sep = "\t")
cancerMF <- as.character(cancerMF$geneset)

CNAmetMF <- read.delim("CNAmet/CNAmet_compPAAD_GOMF/enrichment_results_wg_result1545719558.txt",
                       header = T, sep = "\t")
CNAmetMF <- as.character(CNAmetMF$geneset)

iGCMF <- read.delim("iGC/iGC_compPAAD_GOMF/enrichment_results_wg_result1545721628.txt",
                    header = T, sep = "\t")
iGCMF <- as.character(iGCMF$geneset)

plrsMF <- read.delim("plrs/PLRS_compPAAD_GOMF/enrichment_results_wg_result1545722000.txt",
                     header = T, sep = "\t")
plrsMF <- as.character(plrsMF$geneset)

OncoMF <- read.delim("OncodriveCIS/Oncodrive_compPAAD_GOMF/enrichment_results_wg_result1545721800.txt",
                     header = T, sep = "\t")
OncoMF <- as.character(OncoMF$geneset)

## intersections and Venns
cancerCNAmetMF <- intersect(cancerMF, CNAmetMF)
canceriGCMF <- intersect(cancerMF, iGCMF)
cancerplrsMF <- intersect(cancerMF, plrsMF)
cancerOncoMF <- intersect(cancerMF, OncoMF)

setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns/")
jpeg("PAAD_MFVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetMF,  iGC = iGCMF,  PLRS = plrsMF,
          OncodriveCIS = OncoMF,
          CCGL = cancerMF),
     zcolor = "style", cexil = 1.5, cexsn = 1.5, opacity = 0.3,size = 20)
dev.off()

####################################################################################################

# Reading the files
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/")
cancerPathway <- read.delim("Venns/CancerGenes_Pathway/enrichment_results_wg_result1545723101.txt",
                            header = T, sep = "\t")
cancerPathway <- as.character(cancerPathway$geneset)

CNAmetPathway <- read.delim("CNAmet/CNAmet_compPAAD_Pathway/enrichment_results_wg_result1545720300.txt",
                            header = T, sep = "\t")
CNAmetPathway <- as.character(CNAmetPathway$geneset)

iGCPathway <- read.delim("iGC/iGC_compPAAD_Pathway/enrichment_results_wg_result1545721685.txt",
                         header = T, sep = "\t")
iGCPathway <- as.character(iGCPathway$geneset)

plrsPathway <- read.delim("plrs/PLRS_compPAAD_Pathway/enrichment_results_wg_result1545722012.txt",
                          header = T, sep = "\t")
plrsPathway <- as.character(plrsPathway$geneset)

OncoPathway <- read.delim("OncodriveCIS/Oncodrive_compPAAD_Pathway/enrichment_results_wg_result1545721813.txt",
                          header = T, sep = "\t")
OncoPathway <- as.character(OncoPathway$geneset)

## intersections and Venns
cancerCNAmetPathway <- intersect(cancerPathway, CNAmetPathway)
canceriGCPathway <- intersect(cancerPathway, iGCPathway)
cancerplrsPathway <- intersect(cancerPathway, plrsPathway)
cancerOncoPathway <- intersect(cancerPathway, OncoPathway)

setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/Venns/")
jpeg("PAAD_PathwayVenn.jpg", width = 1260, height = 890) 
venn(list(CNAmet = CNAmetPathway,  iGC = iGCPathway,  PLRS = plrsPathway,
          OncodriveCIS = OncoPathway,
          CCGL = cancerPathway),
     zcolor = "style", cexil = 1.5, cexsn = 1.5, opacity = 0.3,size = 20)
dev.off()
