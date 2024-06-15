#' -----------------------------------------------------------------------------
#' One-sample MR (2SLS; AER)
#' 
#' @description This script performs one-sample MR using individual-level 
#' genotype, methylation and disease status/gene expression data to confirm 
#' directionality of causal estimates for SMR-significant associations 
#' 
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Tues Apr 16 2024. 
#' -----------------------------------------------------------------------------
# load required libraries

library(AER)

# ------------------------------------------------------------------------------
# load individual-level data 

# load dataframe containing individual-level genotype dosages, methylation 
# and disease status/gene expression data

master_df <- readRDS("onesamplemr_dcm_master_df.rds") 
attach(master_df)

# read in CpG-SMR instrumental variable mappings
cpg_IV_map <- readRDS("cpg_IV_map.rds") 

# ------------------------------------------------------------------------------
# loop to perform ivreg for all pairs

regout_ivreg<- matrix(0, nrow(cpg_IV_map), 12)

# nrow(cpg_gene_map)

for (i in 1:nrow(cpg_IV_map)){
  cpg = cpg_IV_map[i,1];
  iv = cpg_IV_map[i,2];
  reg_formula = '';
  reg_formula = as.formula( as.character( paste0("CHF.Etiology", '~', cpg, '|', "`", iv, "`" ) ))
  regout_ivreg[i, 1] <- cpg
  regout_ivreg[i,2] <- iv
  regout_ivreg[i, 3:6] <- summary(ivreg(reg_formula), data= master_df, na.action=na.exclude, diagnostics = T)$coefficient[2,]
  regout_ivreg[i,7:12]<- as.vector(summary(ivreg(reg_formula), data= master_df, na.action=na.exclude, diagnostics = T)$diagnostics[,3:4])
}
colnames(regout_ivreg) <-c("CpG", "smr_IV", "Estimate", "Std.Err", "t-value", "p-value", "weak_instru_stat", "wu_hausman_stat", "sargan_stat",
                           "weak_instru_p", "wu_hausman_p", "sargan_p")

saveRDS(regout_ivreg, "regout_onesamplemr_dcm_smr005.rds")

# ------------------------------------------------------------------------------


