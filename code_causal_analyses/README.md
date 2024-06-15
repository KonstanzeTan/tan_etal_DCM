## General 
This repository contains scripts for primary causal analyses between sentinel CpG methylation and a trait (i.e. DCM/proximal gene expression), one-sample MR validation of SMR-signficant hits, and colocalisation analyses to assess the posterior probability of a shared causal variant underlying the assessed traits. 

## SMR 

The following commands were used to run SMR: 
```
smr-1.3.1 \
    --bfile /1000G_data/1000G_all/1000G_phase3_20130502_chr${chr}_no_duplicates  \ # genotype of all populations in 1000G
    --gwas-summary ./chr_${chr}_gwas_harmonised.ma \	# SNP-DCM associations from GWAS of DCM (UKBB); harmonised with meQTL
    --beqtl-summary ./chr_${chr} \			# BESD-formatted meQTL data  (MAGNet)
    --peqtl-smr 5e-2 \
    --peqtl-heidi 5e-2 \
    --diff-freq 1 \
    --diff-freq-prop 1 \
    --out ./chr_${chr}_DCM >> "${log_file}" 2>&1
```

## One-sample MR

* 'validation_onesampleMR_2SLS.R': performs one-sample MR using individual-level genotype, methylation and disease status/gene expression data to confirm directionality of causal estimates for SMR-significant associations

## Colocalisation analyses 

* 'coloc_meqtl_gwas/eqtl.R': analyses colocalisation between genetic associations of trait pairs (i.e. meQTL and GWAS; meQTL and eQTL). Region for colocalisation
* analysis is selected as +/-500kb of DCM sentinel CpGs, corresponding to the window size used for generating cis-meQTL
