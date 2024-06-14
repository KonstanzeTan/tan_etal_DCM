#' -----------------------------------------------------------------------------
#' Detect consensus modules between ancestry-specific co-methylation networks 
#' (WGCNA)
#' 
#' @description This script reads in ethnic-specific methylation data and 
#' specifies the parameters used to perform consensus network analysis using 
#' WGCNA. It identifies co-methylation moedules that are identified in common 
#' within ancestry-specific analysis 
#' 
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg)
#'
#' @date Thur Mar 7 2024
#' -----------------------------------------------------------------------------
# Load required libraries
suppressMessages({
  library(WGCNA)
  library(cluster)
  library(data.table)
  library(flashClust)
})

# ------------------------------------------------------------------------------

# read in ancestry-specific methylation data for DCM-associated CpGs
## rows = samples, columns = CpGs
datMeth_CAU <- readRDS("datMeth_CAU.rds")
datMeth_AA <- readRDS("datMeth_AA.rds")
# ------------------------------------------------------------------------------
# combine ancestry-specific methylation data for consensus network analysis

# number of networks used in the consensus network analysis 
nSets = 2

# Vector with descriptive names of the two sets
setLabels = c("CAU", "AA")
shortLabels = c("CAU", "AA")

# Define a list whose components contain the data
multiMeth <- vector(mode = "list", length = nSets)
multiMeth[[1]] <- list(data = datMeth_CAU)
names(multiMeth[[1]]$data) <- names(datMeth_CAU)
rownames(multiMeth[[1]]$data) <- dimnames(datMeth_CAU)[[1]]
multiMeth[[2]] <- list(data = datMeth_AA)
names(multiMeth[[2]]$data) <- names(datMeth_AA)
rownames(multiMeth[[2]]$data) <- dimnames(datMeth_AA)[[1]]

# check that the data has the correct format
methSize <- checkSets(multiMeth)

# ------------------------------------------------------------------------------
# run the automatic module detection procedure 
netConsensus <-blockwiseConsensusModules(multiMeth, maxBlockSize = 32198,
                                         power = 12,
                                         minModuleSize = 30,
                                         deepSplit =2, 
                                         pamRespectsDendro = FALSE,
                                         mergeCutHeight = 0.25,
                                         numericLabels = TRUE,
                                         minKMEtoStay = 0,
                                         saveTOMS = TRUE)


# save data
saveRDS(netConsensus, "netConsensus_sft12_sd002_discfdr005_repconstdir.rds")
