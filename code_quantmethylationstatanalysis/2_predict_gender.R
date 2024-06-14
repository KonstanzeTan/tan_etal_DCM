#' -----------------------------------------------------------------------------
#' Check gender of samples based on methylation patterns on sex chromosomes
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Tuesday Aug 22 2023
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# load required libraries

require(minfi)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
require(S4Vectors)

# ------------------------------------------------------------------------------
# read and process methylation array data

setwd("./data/raw_IDAT")
filenames <- unique(substr(dir(), 1, 19))
RGset <- read.metharray(file.path(paste0("./", filenames)), verbose=TRUE)
RGset <- bgcorrect.illumina(RGset)  # Illumina background subtraction 

# ------------------------------------------------------------------------------
# estimate sex of samples 
GMsetEx <- mapToGenome(RGset)
estSex <- getSex(GMsetEx)

# ------------------------------------------------------------------------------
# save results

save(estSex, file="estSex.RData")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()