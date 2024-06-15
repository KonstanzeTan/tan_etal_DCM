#' -----------------------------------------------------------------------------
#' Construction of CVD trait-specific methylation risk score (MRS) in fine-mapped regions
#'
#' @description This script constructs CVD trait-specific MRS for a single
#' fine-mapped region and tests its association with the respective CVD trait. 
#' Weights for MRS construction were obtained from prior association tests between
#' individual CpG loci in a region and CVD traits. 
#' 
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Wed Jun 14 2024
#' -----------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)

# ------------------------------------------------------------------------------
# parse command-line arguments 
args <- commandArgs(trailingOnly = TRUE)

# assign input and output file paths
input_file <- args[1]
output_file_cvd_MRS <- args[2]
output_file_cvd_MRS_assoc <- args[3]

# ------------------------------------------------------------------------------
# load methylation and trait data 
## columns: sample ID, regression covariates, CVD trait information, methylation 
## levels of individual CpGs in region
combined_meth_cvd_covar <- readRDS("combined_meth_cvd_covar.rds")

# ------------------------------------------------------------------------------
# construction of methylation risk score (MRS)

# get weights for MRS construction 
## columns: cpg, beta, se, p, trait
indivcpg_assoc_cvd_500bp <- readRDS("indivcpg_assoc_cvd_500bp.rds")

# subset association results for CpGs in a given region 
ROI_assoc_cvd_500bp <- indivcpg_assoc_cvd_500bp[indivcpg_assoc_cvd_500bp$cpg 
                                                %in% colnames(combined_meth_cvd_covar),]

# ------------------------------------------------------------------------------
# regression of MRS against CVD traits

# specify continuous and binary CVD traits
continuous_cvd_traits <- c("Creat", "DBP", "FramCHD", "SBP", "hsCRP")
binary_cvd_traits <- c("Angina", "CHD", "HT.Hx", "MI", "HT.on.Rx")

# get unique CVD traits 
unique_traits <- unique(ROI_assoc_cvd_500bp$trait)

# initialize dataframe for results
results <- data.frame(trait = character(), beta = numeric(), SE = numeric(), P = numeric(), stringsAsFactors = FALSE)

# loop through each trait to calculate MRS and perform regression
for (trait in unique_traits) {
  # get weights for the current trait (individual loci associations)
  trait_weights <- ROI_assoc_cvd_500bp[ROI_assoc_cvd_500bp$trait == trait, c("cpg", "beta")]
  
  # initialize methylation risk score column
  combined_meth_cvd_data[[paste0("methylation_risk_score_", trait)]] <- 0
  
  # calculate weighted sum for each row
  for (i in 1:nrow(trait_weights)) {
    cpg <- trait_weights$cpg[i]
    weight <- trait_weights$beta[i]
    
    if (cpg %in% colnames(combined_meth_cvd_data)) {
      combined_meth_cvd_data[[paste0("methylation_risk_score_", trait)]] <- 
        combined_meth_cvd_data[[paste0("methylation_risk_score_", trait)]] +
        combined_meth_cvd_data[[cpg]] * weight
    }
  }
  
  # define regression formula
  formula <- as.formula(paste(trait, "~", paste("methylation_risk_score_", trait, " + Age + Sex + CD8T + CD4T + NK + Bcell + Mono + Gran", sep = "")))
  
  # fit model based on trait type
  if (trait %in% continuous_cvd_traits) {
    model <- lm(formula, data = combined_meth_cvd_data)
    summary_model <- summary(model)
    beta <- summary_model$coefficients[paste0("methylation_risk_score_", trait), "Estimate"]
    SE <- summary_model$coefficients[paste0("methylation_risk_score_", trait), "Std. Error"]
    P <- summary_model$coefficients[paste0("methylation_risk_score_", trait), "Pr(>|t|)"]
  } else if (trait %in% binary_cvd_traits) {
    model <- glm(formula, data = combined_meth_cvd_data, family = binomial)
    summary_model <- summary(model)
    beta <- summary_model$coefficients[paste0("methylation_risk_score_", trait), "Estimate"]
    SE <- summary_model$coefficients[paste0("methylation_risk_score_", trait), "Std. Error"]
    P <- summary_model$coefficients[paste0("methylation_risk_score_", trait), "Pr(>|z|)"]
  } else {
    next
  }
  
  # store results to dataframe 
  results <- rbind(results, data.frame(trait = trait, beta = beta, SE = SE, P = P, stringsAsFactors = FALSE))
}

# save results 
saveRDS(combined_meth_cvd_data, output_file_cvd_MRS)
saveRDS(results, output_file_cvd_MRS_assoc)

