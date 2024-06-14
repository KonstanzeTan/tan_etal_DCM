## General 
This repository includes scripts for analyzing expression quantitative trait methylation loci (eQTM) and for estimating hidden sources of variation in gene expression data, which are used as covariates in the analysis.

## cis-eQTM analysis 

* `ciseQTM_analysis.R` : performs cis-eQTM (expression quantitative trait methylation) analysis to investigate
relationships between sentinel CpG methylation and proximal gene expression (<1Mb based on gene transcriptional start site).

## Estimating PEER factors

* 'estimate_peer_factors.R' : learning of hidden determinants in the form of PEER (Probabilistic Estimation of Expression Residuals) from normalized
and preprocessed gene expression data, including model setup and training. 
