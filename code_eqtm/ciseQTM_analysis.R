#' -----------------------------------------------------------------------------
#' cis-eQTM analysis (MatrixeQTL)
#' 
#' @description This script performs cis-eQTM (expression quantitative trait 
#' methylation) analysis to investigate relationships between sentinel CpG methylation
#' and proximal gene expression (<1Mb based on gene transcriptional start site).
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg) 
#'
#' @date Monday Jan 22 2024
#' -----------------------------------------------------------------------------
# load required packages 

library(devtools)
library(Biobase)
library(MatrixEQTL)

# ------------------------------------------------------------------------------
# read in covariate, methylation and gene expression data

cvrt <- read.delim("cvrt_306_AgeGenderRaceRINPEERFACTORS", check.names = FALSE) 
cvrt <- as.data.frame(t(cvrt))  # Row names: covariates (factor encoded), column names: samples
load("beta_306_magnet_sentinels.RData")
gene_expr <- read.delim("magnet_expr_filt_norm", check.names = FALSE)  # Row names: ENSGID, column names: BeadChip

# read in CpG and Gene position information
cpgspos <- read.delim("cpgspos_sentinels") # 
genepos <- read.delim("genepos_16465_auto_hg19")  # Geneid, chr, left, right; autosomal genes with hg19 coordinates

# standardize column order
gene_expr <- gene_expr[, colnames(beta)]
colnames(beta) == colnames(gene_expr)
cvrt <- cvrt[, colnames(beta)]
colnames(beta) == colnames(cvrt)

# prepare SlicedData object for gene expression
gene <- SlicedData$new()
gene$CreateFromMatrix(as.matrix(gene_expr))
gene$fileSliceSize <- 2000  # Read file in pieces of 2,000 rows

# prepare SlicedData object for covariates
cvrt <- SlicedData$new()  # Leave as this if not running with covariates 
cvrt$CreateFromMatrix(as.matrix(cvrt))

# prepare SlicedData object for CpGs (matrix, rownames=CpGs, colnames=samples)
cpgs <- SlicedData$new()
cpgs$CreateFromMatrix(as.matrix(beta))
cpgs$fileSliceSize <- 2000  # Read file in pieces of 2,000 rows

output_file_name <- "eQTM_AgeGenderRaceRINPEERFACTORS_cis_1Mb"

me <- Matrix_eQTL_main(
  snps = cpgs,  # Sliced data object 
  gene = gene,
  cvrt = cvrt,
  output_file_name.cis = output_file_name,
  pvOutputThreshold.cis = 1,
  pvOutputThreshold = 0,  # Set = 0 if analyzing cis-pairs
  snpspos = cpgspos,
  genepos = genepos,
  useModel = modelLINEAR,  # standard additive linear model rather than dominant or recessive linear model (applies to genetic data)
  errorCovariance = numeric(),  # Independence for gene expression
  verbose = TRUE,
  cisDist = 1e6,  # Distance for local gene-CpG pairs. 
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE  # Calculate FDR?
)
