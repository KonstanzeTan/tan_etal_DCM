#' -----------------------------------------------------------------------------
#' Check for duplicate samples using methylation betas at SNPs
#' 
#' @description This script retrieves methylation betas at 59 SNPs on the EPIC 
#' array. Duplicate samples are identified as those with high correlation (
#' r>0.8, Pearson correlation) 
#' in methylation betas across the SNPs. 
#' 
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Tuesday July 11 2023
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# load required libraries

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(tidyverse)
library(Hmisc)
library(gplots)

# ------------------------------------------------------------------------------
# Retrieve methylation betas at 59 SNPs 

# set working directory to IDAT files directory
setwd("./data/raw_IDAT")

# read IDAT files 
filenames <- unique(substr(dir(),1,19))

RGset <- read.metharray(file.path(paste0("./", filenames)), verbose=TRUE)
RGset <- bgcorrect.illumina(RGset)

snps <- getSnpBeta(RGset)
head(snps)

save(snps, file="snps_beta.RData")
# ------------------------------------------------------------------------------
# Check for duplicate samples 

# calculate correlation coefficients and p-values 
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# perform duplicates check
res <- rcorr(snps, type="pearson")
res_flat <- flattenCorrMatrix(res$r, res$P)
res_high <- res_flat[res_flat$cor > 0.8, ]

# save table of highly correlated samples 
write.table(res_high, "highly_correlated_samples.txt", quote=F, sep = "\t", 
            row.names=F, col.names=T)
# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()