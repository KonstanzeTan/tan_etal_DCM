#' -----------------------------------------------------------------------------
#' Conditional regression to identify independent signals 
#' 
#' @description This script was used to identify independent signals amongst 
#' candidate sentinel CpGs at genomic loci with >1 CpG associated with DCM 
#' (discovery EWAS Bonferroni-corrected P<0.05). Signals were considered 
#' sentinel CpGs if they remained significantly associated with DCM 
#' (Bonferroni-corrected P<0.05) after conditioning on the lead signal. As in EWAS,
#' conditional analysis was first conducted separately by ancestry, followed by 
#' trans-ancestry meta-analysis. Conditional analysis in the African-American (AA)
#' ancestry group is used as an example. 
#' 
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg>
#'
#' @date Wed Mar 6 2024
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Create ancestry-specific phenotype files 

phe<-read.delim("phe_genderpass_rmDup.txt") 

## split phe by Ancestry
phe_AA<-phe[phe$Race=="AA",]
dim(phe_AA) 
phe_AA$DCM_status <- factor(phe_AA$DCM_status, levels = c("NF", "ID-CMP")) # ID-CMP = idiopathic DCM; NF = non-failing
phe_AA$Gender <- factor(phe_AA$Gender, levels = c("Male", "Female"))

# ------------------------------------------------------------------------------
# load race-specific control-probe PCs
load("PC_cp_AA.RData")

# ------------------------------------------------------------------------------
# read in conditional analysis CpG pairs

cond_analysis_pairs<-readRDS("cond_analysis_pairs.rds")
cpg_id<-unique(c(cond_analysis_pairs[,2], cond_analysis_pairs[,3]))

# ------------------------------------------------------------------------------
# Get ancestry-specific quantile-normalised methylation betas 

##  load quantile-normalised betas 
beta <- readRDS("beta_QN_detP001_marker095.rds")
dim(beta) 

## get ancestry-specific betas for candidate sentinel CpGs 
beta_AA <- beta[rownames(beta) %in% cpg_id,colnames(beta) %in% phe_AA$BeadChip]
beta_AA<-as.data.frame(t(beta_AA))
beta_AA$BeadChip<-rownames(beta_AA)

# ------------------------------------------------------------------------------                       
# Create master dataframe of covariates for conditional analysis

## format control probe PCs
PC_cp_AA<-as.data.frame(PC_cp_AA)
PC_cp_AA$BeadChip<- row.names(PC_cp_AA)
regout_AA_conditional_all_variables<-merge(phe_AA, PC_cp_AA, by = "BeadChip")
regout_AA_conditional_all_variables<-merge(regout_AA_conditional_all_variables, beta_AA, by="BeadChip", all.x=TRUE)

## create independent variable list (covariates)
Header_List= c(); 
Header_List= colnames(regout_AA_conditional_all_variables)
PC_cp_List = c(Header_List[15:24])
print(PC_cp_List)
RHS_equation = 'Age+ as.factor(Gender)'
for (element_loop in PC_cp_List) {
  
  RHS_equation = paste(RHS_equation, '+', element_loop, sep="")
  
}
print(RHS_equation)

## check factors levels
levels(as.factor(regout_AA_conditional_all_variables$Gender))
levels(regout_AA_conditional_all_variables$DCM_status)

# ------------------------------------------------------------------------------
# Conditional analysis

regout_AA_conditional<- matrix(0, nrow(cond_analysis_pairs), 4)

for (i in 1:nrow(cond_analysis_pairs)){
  lead_cpg = cond_analysis_pairs[i,2];
  nonlead_cpg = cond_analysis_pairs[i,3];
  reg_formula = '';
  reg_formula = as.formula( as.character( paste('DCM_status', '~', nonlead_cpg,
                                                '+', RHS_equation,'+', lead_cpg ) ))
  regout_AA_conditional[i,] <- summary(glm(reg_formula, family= binomial, 
                                           data= regout_AA_conditional_all_variables, 
                                           na.action=na.exclude))$coefficient[2,]}

colnames(regout_AA_conditional)<-c("Estimate", "Std_Error", "t_value", "P")
regout_AA_conditional<- as.data.frame(regout_AA_conditional)
regout_AA_conditional$loci_ID <- cond_analysis_pairs$loci_ID
regout_AA_conditional$lead_cpg <- cond_analysis_pairs$lead_cpg
regout_AA_conditional$nonlead_cpg <- cond_analysis_pairs$nonlead_cpg

saveRDS(regout_AA_conditional, "regout_AA_conditional.rds")

# ------------------------------------------------------------------------------
sessionInfo()
