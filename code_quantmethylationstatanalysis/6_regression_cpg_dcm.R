#' -----------------------------------------------------------------------------
#' Association testing between individual CpG loci and DCM
#' 
#' @description This script was used for discovery 
#' and replication-stage association testing between individual
#' CpG loci and DCM. Discovery-stage EWAS (MAGNet) was performed separately by 
#' ancestry (African American, Caucasian). Ancestry-specific EWAS of DCM using
#' methylation data from African Americans (AA) is used as an example. 
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
# Calculate ancestry-specific control probe principal components (PCs) 

load("ctrlprobes.RData")

## split ctrlprobes by ancestry
ctrl.all_AA <- ctrl.all[rownames(ctrl.all) %in% phe_AA$BeadChip,]

## calculate principal components
pca_cp_AA <- prcomp(na.omit(ctrl.all_AA))
PC_cp_AA <- predict(pca_cp_AA)
dim(PC_cp_AA)
save(PC_cp_AA, file="PC_cp_AA.RData")

# ------------------------------------------------------------------------------
# Get ancestry-specific quantile-normalised methylation betas 

##  load quantile-normalised betas 
beta <- readRDS("beta_QN_detP001_marker095.rds")
dim(beta) 

## select ancestry-specific betas
beta_AA <- beta[,colnames(beta) %in% phe_AA$BeadChip]
dim(beta_AA) 
# check that there are no rows with all NA in the ethnic-specific betas
num_rows_with_na <- sum(rowSums(is.na(beta_AA)) == ncol(beta_AA)

# ------------------------------------------------------------------------------
#  Create master dataframe of covariates for EWAS 

## format control probe PCs
PC_cp_AA<-as.data.frame(PC_cp_AA)
PC_cp_AA$BeadChip<- row.names(PC_cp_AA)
phe_AA_all_variables<-merge(phe_AA, PC_cp_AA, by = "BeadChip")
attach(phe_AA_all_variables)

## ensure matching order of samples between covariates and methylation beta df
beta_order<-phe_AA_all_variables$BeadChip
beta_AA <- beta_AA[,match(beta_order, colnames(beta_AA))]
colnames(beta_AA)==phe_AA_all_variables$BeadChip

## create independent variable list (covariates)
Header_List= c(); 
Header_List= colnames(phe_AA_all_variables)
Independent_Variable_List = c(Header_List[15:24])
print(Independent_Variable_List)
RHS_equation = 'beta_AA[i,] +Age+ as.factor(Gender)'
for (element_loop in Independent_Variable_List) {
  
  RHS_equation = paste(RHS_equation, '+', element_loop, sep="")
  
}
print(RHS_equation)

## check factor levels
levels(as.factor(Gender))
levels(DCM_status)

# ------------------------------------------------------------------------------
# Logistic regression of each CpG site against DCM

regout_AA<- matrix(0, nrow(beta_AA), 4)

for (i in 1:nrow(beta_AA)){
  reg_formula = '';
  reg_formula = as.formula( as.character( paste('DCM_status~', RHS_equation, sep="") ) )
  regout_AA[i,] <- summary(glm(reg_formula,family=binomial, na.action=na.exclude))$coefficient[2,]}

rownames(regout_AA) <- rownames(beta_AA)
dim(regout_AA)
print(reg_formula)
saveRDS(regout_AA, file="regout_AA_AgeGender10PCs_marker095.rds")

# ------------------------------------------------------------------------------
sessionInfo()