MethylMix_ModelGeneExpression_FDR <- function (METcancer, GEcancer, CovariateData = NULL) 
{
  OverlapSamples = intersect(colnames(METcancer), colnames(GEcancer))
  cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
  GEcancer = GEcancer[, OverlapSamples, drop = FALSE]
  METcancer = METcancer[, OverlapSamples, drop = FALSE]
  if (!is.null(CovariateData)) 
    CovariateData = CovariateData[OverlapSamples, , drop = FALSE]
  Rsquares = matrix(0, nrow = length(rownames(METcancer)), 
                    ncol = 1)
  Genes = rownames(METcancer)
  PvalueThreshold = 0.001
  RsquareThreshold = 0.1
  cat("Correlating methylation data with gene expression...\n")
  i <- NULL
  Rsquares = foreach::foreach(i = 1:length(rownames(METcancer)), 
                              .combine = "c") %dopar% {
                                Rsq = c(0,0)
                                tmpGene = unlist(strsplit(Genes[i], "---"))[1]
                                pos = which(rownames(GEcancer) == tmpGene)
                                if (length(pos) > 0) {
                                  if (!is.null(CovariateData)) {
                                    res = lm(GEcancer[pos, ] ~ METcancer[Genes[i], ] + factor(CovariateData))
                                    res.summary = summary(res)
                                    an = anova(res)
                                    if (res$coefficients[2] < 0 & res.summary$coefficients[2,4] < PvalueThreshold) {
                                      Rsq = res.summary$r.squared
                                    }
                                  }
                                  else {
                                    res = lm(GEcancer[pos, ] ~ METcancer[Genes[i], ])
                                    res.summary = summary(res)
                                    if (res$coefficients[2] < 0 & res.summary$coefficients[2,4] < PvalueThreshold) {
                                      Rsq = c(res.summary$r.squared, res.summary$coefficients[2,4])
                                    }
                                  }
                                }
                                Rsq
                              }
  Rsquares = data.frame(matrix(Rsquares, ncol = 2, byrow = T))
  Rsquares[,3] = p.adjust(Rsquares[,2], method = "fdr", n = length(Rsquares[,2])) # estimating fdr
  FunctionalGenes = Genes[Rsquares[,1] > RsquareThreshold]
  cat("\nFound", length(FunctionalGenes), "transcriptionally predictive genes.\n")
  
  Rsquares = Rsquares[Rsquares[,1] > RsquareThreshold,]
  FunctionalGenes_withfdr = cbind(FunctionalGenes, Rsquares[,3])
  return(FunctionalGenes_withfdr)
}

