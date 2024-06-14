#' -----------------------------------------------------------------------------
#' Estimate PEER factors from gene expression data
#' 
#' @description: This script learns hidden determinants in the form of PEER 
#' (Probabilistic Estimation of Expression Residuals) from normalized
#' and preprocessed gene expression data, including model setup and training. 
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Mon Jan 29 2024
#' -----------------------------------------------------------------------------
# load required libraries 

library(peer)
library(ggplot2)

# ------------------------------------------------------------------------------
# read in expression and covariate data 

# read normalised, preprocessed expression data 
expr <- read.delim("magnet_expr_filt_norm", check.names = FALSE)
dim(expr)

# Format for PEER factor learning
expr <- t(expr)  # Rows = samples, columns = genes

# read in covariates file
covs <- read.delim("cvrt_306_dcmstatusAgeGenderRaceRIN")
covs[] <- lapply(covs, as.numeric)
covs <- covs[rownames(expr), ]
stopifnot(rownames(covs) == rownames(expr))

# ------------------------------------------------------------------------------
# Initialize and set up the PEER model

model <- PEER()
PEER_setPhenoMean(model, as.matrix(expr))
dim(PEER_getPhenoMean(model))

# Set number of PEER factors
PEER_setNk(model, 35)  
PEER_getNk(model)

# Additional PEER model settings
PEER_setAdd_mean(model, TRUE)
PEER_setNmax_iterations(model, 10000)
PEER_setCovariates(model, as.matrix(covs))

# Train the PEER model
PEER_update(model)

# Extract results from the PEER model
PEER_factors <- PEER_getX(model)  # Posterior mean of the inferred confounders (NxK matrix)
dim(PEER_factors)
PEER_weights <- PEER_getW(model)
PEER_bounds <- PEER_getBounds(model)
PEER_vars <- PEER_getResidualVars(model)
PEER_alpha <- PEER_getAlpha(model)

rownames(PEER_factors) <- rownames(expr)

# ------------------------------------------------------------------------------
# Write results to files

output_prefix <- "PEER_magnet_dcmAgeGenderRaceRIN"
write.table(PEER_factors, paste0(output_prefix, "_factors"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(PEER_bounds, paste0(output_prefix, "_bounds"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(PEER_vars, paste0(output_prefix, "_vars"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(PEER_alpha, paste0(output_prefix, "_alpha"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(PEER_weights, paste0(output_prefix, "_weights"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# ------------------------------------------------------------------------------
# Plot factor relevance to decide on number of PEER factors to include 
# (PEER factors before the 'elbow' inflection point are included)

PEER_alpha$variance <- 1/PEER_alpha$V1
PEER_alpha$index <- 1:nrow(PEER_alpha)
PEER_alpha <- PEER_alpha[order(PEER_alpha$variance, decreasing=T),]
 
ggplot(PEER_alpha, aes(x = seq_along(variance), y = variance)) +
  geom_line(color = "blue") +  # Connect the points with lines
  geom_point(color = "red", size = 2) +  # Add red points
  geom_text(aes(label = index), vjust = -0.5, hjust = 1) +  # Add text labels
  labs(title = "Connected Scatterplot for Variance",
       x = "Observation",
       y = "Variance")

sessionInfo()
