#' -----------------------------------------------------------------------------
#' Colocalisation analysis of genetic assoiciations (coloc)
#' 
#' @description This script analyses colocalisation between genetic associations 
#' of trait pairs (i.e. meQTL and GWAS; meQTL and eQTL). Region for colocalisation
#' analysis is selected as +/-500kb of DCM sentinel CpGs, corresponding to the 
#' window size used for generating cis-meQTL
#' 
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg> 
#'
#' @date Tues Apr 16 2024
#' -----------------------------------------------------------------------------
# load required libraries

library(coloc)
# ------------------------------------------------------------------------------
# read in gwas and meqtl data
## large list containing invidual list of cpg-specific meqtl/gwas associations
## each individual list contains beta, varbeta, snp, position, CpG and type
## e.g. of list names (cg12359658_meqtl, cg12359658_gwas)

coloc_meqtl_gwas_data <- readRDS("coloc_meqtl_gwas_data.rds") 
attach(coloc_meqtl_gwas_data)

# read in coloc file assignments: pairs of meqtl and gwas datasets to be analysed
## dataframe containing: cpg, meqtl_file and gwas_file columns 
## each column contains the names of lists in coloc_meqtl_gwas_data
## (e.g. cg12359658_meqtl, cg12359658_gwas)

coloc_assignments<- readRDS("coloc_assignments.rds")

# ------------------------------------------------------------------------------
# run coloc analysis for all cpg-gene pairs with smr-sig results: 

# Initialize the master dataframe
results_coloc<- matrix(0, nrow(coloc_assignments), 7)

# Loop through each row in smr_sig_cpg_gene
for (i in 1:nrow(coloc_assignments)) {
  cpg <- coloc_assignments$cpg[i]
  meqtl_file <-  coloc_assignments$meqtl_file[i]
  gwas_file <- coloc_assignments$gwas_file[i]
  meqtl_dataset <-coloc_meqtl_gwas_data[[meqtl_file]]
  gwas_dataset <- coloc_meqtl_gwas_data[[gwas_file]]
  results_coloc[i,2:7] <- coloc.abf(meqtl_dataset, gwas_dataset)$summary
  results_coloc[i,1]<- cpg
}

results_coloc <- as.data.frame(results_coloc)
names(results_coloc) <- c("cpg", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

# save results
saveRDS(results_coloc, "results_coloc.rds")


