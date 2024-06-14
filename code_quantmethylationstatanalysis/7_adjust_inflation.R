#' -----------------------------------------------------------------------------
#' Estimating and adjusting for test statistic inflation (BACON)
#' 
#' @description A Bayesian method implemented within the BACON R package was
#' applied to estimate test statistic inflation in ancestry-specific EWAS 
#' and trans-ancestry meta-analysis. This script is run to
#' estimate inflation and adjust effect sizes and standard errors accordingly
#' EWAS performed in the African-American ancestry group is used as an example.
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Wed Mar 6 2024
#' -----------------------------------------------------------------------------
#' 
# ------------------------------------------------------------------------------
# load required libraries 

library(bacon)

# ------------------------------------------------------------------------------
# format data for BACON

# load EWAS results
regout_AA <- readRDS("regout_AA_AgeGender10PCs_marker095.rds")
colnames(regout_AA)<-c('Estimate','Std. Error', 'z value','Pr(>|z|)')

# create separate matrices for effect size and standard error
AA_beta<-as.matrix(as.numeric(regout_AA[,1]))
rownames(AA_beta)<-rownames(regout_AA)
AA_se<-as.matrix(as.numeric(regout_AA[,2]))
rownames(AA_se)<-rownames(regout_AA)
bc<-bacon(NULL, AA_beta, AA_se)
bc

# plot Gibbs Sample mixture fit 
fit(bc, n = 100)
# ------------------------------------------------------------------------------
# Get inflation-corrected effect sizes and SE for trans-ancestry meta-analysis

Beta_baconcorr<-as.data.frame(es(bc, corrected=TRUE))
Beta_baconcorr$CpG<-row.names(Beta_baconcorr)
names(Beta_baconcorr)<-c("Beta_baconcorr", "CpG")
SE_baconcorr<-as.data.frame(se(bc, corrected=TRUE))
SE_baconcorr$CpG<-row.names(SE_baconcorr)
names(SE_baconcorr)<-c("SE_baconcorr", "CpG")
regout_AA<-as.data.frame(regout_AA)
regout_AA$CpG<-row.names(regout_AA)
regout_AA<-merge(regout_AA, Beta_baconcorr, by="CpG", all.x=TRUE)
regout_AA<-merge(regout_AA, SE_baconcorr, by="CpG", all.x=TRUE)

# save file 
saveRDS(regout_CAU, "regout_CAU_baconcorr.rds")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()