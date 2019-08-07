HUGOtoEntrez <- function(InputGenes){
  PCDat= PCAliasDat = PCPreviousDat = data.frame(matrix(nrow = 0, ncol = 2))
  NCDat= NCAliasDat = NCPreviousDat = data.frame(matrix(nrow = 0, ncol = 2))
  PhenoDat= PhenoAliasDat = PhenoPreviousDat = data.frame(matrix(nrow = 0, ncol = 2))
  PseudoDat= PseudoAliasDat = PseudoPreviousDat = data.frame(matrix(nrow = 0, ncol = 2))
  OtherDat= OtherAliasDat = OtherPreviousDat = data.frame(matrix(nrow = 0, ncol = 2))
  
  ApprovedEntrez <- function(Genes, Data){
    Dat <- data.frame(matrix(ncol = 2, nrow = length(Genes)))
    colnames(Dat) <- c("HUGO", "Entrez")
    if(length(Genes)!=0){
      for (i in 1:length(Genes)) {
        Dat[i,1] <- Genes[i]
        Dat[i,2] <- Data[which(Data$Approved.symbol==Genes[i]), 5][1]
      }
    }
    return(Dat)
  }
  AliasEntrez <- function(Genes, Data){
    Dat <- data.frame(matrix(ncol = 2, nrow = length(Genes)))
    colnames(Dat) <- c("HUGO", "Entrez")
    if(length(Genes)!=0){
      for (i in 1:length(Genes)) {
        Dat[i,1] <- Genes[i]
        Dat[i,2] <- Data[which(Data$Alias.symbol==Genes[i]), 5][1]
      }
    }
    return(Dat)
  }
  PreviousEntrez <- function(Genes, Data){
    Dat <- data.frame(matrix(ncol = 2, nrow = length(Genes)))
    colnames(Dat) <- c("HUGO", "Entrez")
    if(length(Genes)!=0){
      for (i in 1:length(Genes)) {
        Dat[i,1] <- Genes[i]
        Dat[i,2] <- Data[which(Data$Previous.symbol==Genes[i]), 5][1]
      }
    }
    return(Dat)
  }
  
  setwd("/home/anita/Benchmarking/two_omics/")
  pc <- read.delim("ProteinCodingGenesHUGO.txt", header = T, sep = "\t")
  nc <- read.delim("NonCodingGenesHUGO.txt", header = T, sep = "\t")
  pheno <- read.delim("PhenotypeGenesHUGO.txt", header = T, sep = "\t")
  pseudo <- read.delim("PseudoGenesHUGO.txt", header = T, sep = "\t")
  other <- read.delim("OtherGenesHUGO.txt", header = , sep = "\t")
  
  allApproved <- unique(c(as.character(pc$Approved.symbol), as.character(nc$Approved.symbol),
                  as.character(pheno$Approved.symbol), as.character(pseudo$Approved.symbol),
                  as.character(other$Approved.symbol)))
  
  allAlias <- unique(c(as.character(pc$Alias.symbol), as.character(nc$Alias.symbol),
                          as.character(pheno$Alias.symbol), as.character(pseudo$Alias.symbol),
                          as.character(other$Alias.symbol)))
  
  allPrevious <- unique(c(as.character(pc$Previous.symbol), as.character(nc$Previous.symbol),
                          as.character(pheno$Previous.symbol), as.character(pseudo$Previous.symbol),
                          as.character(other$Previous.symbol)))
  
  ApprovedGenes <- intersect(InputGenes, allApproved)
  if (length(ApprovedGenes)!=0){
    a <- intersect(ApprovedGenes, as.character(pc$Approved.symbol))
    b <- intersect(ApprovedGenes, as.character(nc$Approved.symbol))
    c <- intersect(ApprovedGenes, as.character(pheno$Approved.symbol))
    d <- intersect(ApprovedGenes, as.character(pseudo$Approved.symbol))
    e <- intersect(ApprovedGenes, as.character(other$Approved.symbol))
    
    PCDat <- ApprovedEntrez(Genes = a, Data = pc)
    NCDat <- ApprovedEntrez(Genes = b, Data = nc)
    PhenoDat <- ApprovedEntrez(Genes = c, Data = pheno)
    PseudoDat <- ApprovedEntrez(Genes = d, Data = pseudo)
    OtherDat <- ApprovedEntrez(Genes = e, Data = other)
  }
  
  RemainGenes <- setdiff(InputGenes, ApprovedGenes)
  AliasGenes <- intersect(RemainGenes, allAlias)
  if(length(AliasGenes)!=0){
    a <- intersect(AliasGenes, as.character(pc$Alias.symbol))
    b <- intersect(AliasGenes, as.character(nc$Alias.symbol))
    c <- intersect(AliasGenes, as.character(pheno$Alias.symbol))
    d <- intersect(AliasGenes, as.character(pseudo$Alias.symbol))
    e <- intersect(AliasGenes, as.character(other$Alias.symbol))
    
    PCAliasDat <- AliasEntrez(Genes = a, Data = pc)
    NCAliasDat <- AliasEntrez(Genes = b, Data = nc)
    PhenoAliasDat <- AliasEntrez(Genes = c, Data = pheno)
    PseudoAliasDat <- AliasEntrez(Genes = d, Data = pseudo)
    OtherAliasDat <- AliasEntrez(Genes = e, Data = other)
  }
  
  RemainGenes2 <- setdiff(InputGenes, c(ApprovedGenes, AliasGenes))
  PreviousGenes <- intersect(RemainGenes2, allPrevious)
  if(length(PreviousGenes)!=0){
    a <- intersect(PreviousGenes, as.character(pc$Previous.symbol))
    b <- intersect(PreviousGenes, as.character(nc$Previous.symbol))
    c <- intersect(PreviousGenes, as.character(pheno$Previous.symbol))
    d <- intersect(PreviousGenes, as.character(pseudo$Previous.symbol))
    e <- intersect(PreviousGenes, as.character(other$Previous.symbol))
    
    PCPreviousDat <- PreviousEntrez(Genes = a, Data = pc)
    NCPreviousDat <- PreviousEntrez(Genes = b, Data = nc)
    PhenoPreviousDat <- PreviousEntrez(Genes = c, Data = pheno)
    PseudoPreviousDat <- PreviousEntrez(Genes = d, Data = pseudo)
    OtherPreviousDat <- PreviousEntrez(Genes = e, Data = other)
  }
  
  
  return(rbind(PCDat, PCAliasDat, PCPreviousDat, NCDat, NCAliasDat, NCPreviousDat,
               PhenoDat, PhenoAliasDat, PhenoPreviousDat, PseudoDat, PseudoAliasDat, PseudoPreviousDat,
               OtherDat, OtherAliasDat, OtherPreviousDat))
}
