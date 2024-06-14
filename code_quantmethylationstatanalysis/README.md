## General 
The scripts in this various repository facilitate tasks in methylation data preprocessing and quality control, calculation of methylation betas, epigenome-wide association analysis, meta-analysis, adjustment for test-statistic inflation and finally the identification of sentinel CpGs. 

## Methylation data preprocessing

* `1_get_intensities.R` : retrieves and normalizes methylation data from Illumina Human Methylation EPIC arrays. It includes background correction, separation and extraction of Type I and Type II probe intensities, and calculation of detection p-values. Additionally, it performs PCA on control probe intensities.
* `2_predict_gender.R` : predicts gender of samples based on methylation patterns on sex chromosomes.
* `3_check_gender_swap.R` : flags samples with mismatched self-reported and predicted gender for exclusion
* `4_check_duplicates.R` : retrieves methylation betas at 59 SNPs on the EPIC array. Duplicate samples are identified as those with high correlation (r>0.8, Pearson correlation) in methylation betas across the SNPs

## Calculating methylation betas
* `5_get_methylation_betas.R` : performs quantile normalisation of marker intensity values and uses the normalised values to calculate the beta value (methylation level at each CpG site). Additionally, a call rate threshold is applied to select samples and markers for further analysis

## Epigenome-wide association analysis of DCM 
* `6_regression_cpg_dcm.R` : perform association testing between individual CpG loci and DCM. Discovery-stage EWAS (MAGNet) was performed separately by ancestry (African American, Caucasian). Ancestry-specific EWAS of DCM using methylation data from African Americans (AA) is used as an example in this script.
* `7_adjust_inflation.R` : utilises a Bayesian method implemented in the bacon R package to estimate test statistic inflation in ancestry-specific EWAS and trans-ancestry meta-analysis and to adjust effect sizes and standard errors accordingly.

## Trans-ancestry meta-analysis (METAL)

* `8_metaanalyse_transancestry_metal.txt` : this is the runfile for inverse-variance weighted meta-analysis of ancestry-specific regression results of individual CpG loci against DCM (CAU, AA). 

To run this script, set the working directory to the same directory as where the METAL software resides. Next, 
execute the runfile as follows 
```
SOURCE 8_metaanalyse_transancestry_metal.txt
```

## Identification of sentinel CpGs

* `9_identify_sentinelcpgs_dcm.R` : assigns candidate sentinel CpGs to unique genomic loci (>1Mb apart).
* `10_regression_conditional.R` : identifies independent signals amongst candidate sentinel CpGs at genomic loci with >1 CpG associated with DCM (discovery EWAS Bonferroni-corrected P<0.05). Signals were considered sentinel CpGs if they remained significantly associated with DCM (Bonferroni-corrected P<0.05) after conditioning on the lead signal. As in EWAS, conditional analysis was first conducted separately by ancestry, followed by trans-ancestry meta-analysis. 
