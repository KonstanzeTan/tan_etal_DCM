#' -----------------------------------------------------------------------------
#' Identify samples with mismatch between predicted and reported gender
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Tuesday Jul 11 2023
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# load required data

# predicted gender 
load("estSex.RData")

# phenotype file containing reported gender 
phe <- read.delim("phe.txt")
# ------------------------------------------------------------------------------
# compare predicted against reported gender 

estSex <- as.data.frame(estSex)
estSex<-estSex[rownames(estSex) %in% phe$BeadChip,]
estSex$BeadChip<-row.names(estSex)
estSex<-merge(estSex, phe[,c("BeadChip", "Gender")], by="BeadChip", all.x=TRUE)

# relabel gender as "F" and "M" for comparison
estSex$Gender[estSex$Gender == "Female"] <- "F"
estSex$Gender[estSex$Gender == "Male"] <- "M"
estSex$Gender_predictedSex_Match <- ifelse(estSex$Gender == estSex$predictedSex, "pass", "fail")

# ------------------------------------------------------------------------------
# save results

write.table(estSex, file="gender_check_status.txt",
            row.names=F, col.names = T, sep="\t", quote=F)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()