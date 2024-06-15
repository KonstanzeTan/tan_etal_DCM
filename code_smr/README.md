## General 
This repository contains scripts for primary causal analysis using SMR, validation using one-sample MR, and colocalisation testing to assess the
posterior probability of a single causal variant underlying the assessed traits


## SMR

SMR was performed on chromosome-specific genetic association and genotype data. Following which, results from individual chromosomes were combined
```
smr-1.3.1 \
    --bfile /1000G_data/1000G_all/1000G_phase3_20130502_chr${chr}_no_duplicates  \ # genotype of all populations in 1000G
    --gwas-summary ./chr_${chr}_gwas_harmonised.ma \	# SNP-DCM associations from GWAS of DCM (UKBB); harmonised with meQTL
    --beqtl-summary ./chr_${chr} \			# BESD-fomratted meQTL data  (MAGNet)
    --peqtl-smr 5e-2 \
    --peqtl-heidi 5e-2 \
    --diff-freq 1 \
    --diff-freq-prop 1 \
    --out ./chr_${chr}_DCM >> "${log_file}" 2>&1

```

## One-sample MR 

One-sample MR was performed on SMR-significant associations to confirm the directionality of causal estimates for SMR-significant associations 

* 'validation_onesampleMR_2SLS.R': performs one-sample MR using individual-level genotype, methylation and disease status/gene expression data 

## Colocalisation 

Colocalisation of meQTl (MAGNet) and GWAS/eQTL (public datasets) were examined for SMR-significant pairs.

* 'coloc_meqtl_gwas/eqtl.R': analyses colocalisation between genetic associations of trait pairs (i.e. meQTL and GWAS; meQTL and eQTL). Region for colocalisation
analysis is selected as +/-500kb of DCM sentinel CpGs, corresponding to the window size used for generating cis-meQTL



