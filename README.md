# Code for Epigenome-wide association analysis Tan et el. 2024

Citation: 
**[title.](URL)**

Authors

## Table of Contents


   * [Identification of meQTLs](#identification-of-meqtls)
      * [EPIC meQTL](#epic-meqtl)
      * [Conditional analysis and pruning](#conditional-analysis-and-pruning)
   * [meQTL replication meDIPseq](#meqtl-replication-medipseq) 
   * [Identification of eQTLs and eQTMs](#identification-of-eqtls-and-eqtms)
   * [Functional analyses of meQTL CpGs](#functional-analyses-of-meqtl-cpgs)
   * [Functional analyses of meQTL pairs](#functional-analyses-of-meqtl-pairs)
      * [Chromatin state enrichment](#chromatin-state-enrichment)
      * [HiC enrichment](#hic-enrichment)
      * [Regulator enrichment](#regulator-enrichment)
      * [Enrichment of transcription factor binding sites at trans-meQTL CpG sites](#enrichment-of-transcription-factor-binding-sites-at-trans-meqtl-cpg-sites)
      * [Network analysis](#network-analysis)
   * [Identification of iQTL](#identification-of-iqtl)
   * [EPIC comparison](#epic-comparison)
   * [Docker](#docker)


## Quantification and statistical analysis of DNA methylation 

We ran EWAS of DCM in two independent cohorts to identify conditionally-independent sentinel CpGs (P<5.96E-08) with confirmed consistent directionality of association. 

The scripts in this repository facilitate various tasks in (i) methylation data preprocessing: including extracting raw intensities, predicting gender, checking for errors and duplicates, calculating methylation beta values, and (ii) Epigenome-wide association analysis: including performing ancestry-specific regression and trans-ancestry meta-analysis, adjusting for test statistic inflation, and the identification of sentinel CpGs marking distinct genomic loci.

The model used in EWAS was: 

```
DCM ~beta(quantile-normalised) + Age + Gender + principal components of control probes capturing ~95% of control probe variation 
```

Software/packages used: minfi (R), bacon (R)

## Enrichment analysis of genomic regulatory features

We analysed sentinel CpGs for enrichment in various genomic regulatory features using a permutation testing approach. 

The scripts in this repository outline the manual construction of a background for enrichment analysis, consisting of 1,000 sets of EPIC array CpGs matched to sentinel CpGs by methylation levels and variability. Additionally, they detail the calculation of permutation p-values to assess enrichment.

Software/packages used: local installation of eFORGE v2 downloaded from [Altius Institute eForge](https://eforge.altiusinstitute.org/)

## Expression quantitative trait methylation analysis

We identified significant associations between sentinel CpG methylation and proximal gene expression (+/- 1Mb) of the sentinel CpG. Replicated eQTMs were defined as significant 
eQTMs in the discovery cohort (MAGNet; disc FDR P<0.05) for which replication testing (BMCB) confirmed consistent directionality of association. 
```
gene expression ~  methylation + Age + Gender + ancestry + RNA Integrity Number(RIN) + 5 PEER factors
```
Software/packages used: MatrixEQTL (R)

## Causal analyses (Mendelian Randomisation and Colocalisation) 

We elucidated the putative causal contribution of sentinel CpG methylation to DCM and proximal gene expression using [SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview).

The scripts in this repository detail the primary SMR analysis, one-sample MR validation of SMR-significant hits, and colocalisation analysis of genetic associations for the assessed traits

Software/packages used: SMR, AER(R), coloc(R)

## Weighted gene correlation network analysis (WGCNA)

We employed WGCNA to detect co-methylation modules amongst DCM-associated CpGs.

Software/packages used: WGCNA (R)

## Construction of methylation risk score (MRS) from fine-mapping data and examining associations with CVD

Following fine-mapping of 28 DCM sentinel CpGs and testing individual CpGs within fine-mapped regions for CVD trait associations, we constructed CVD-trait specific MRS as a weighted sum of methylation betas from CpGs within fine-mapped regions and tested MRS for association with CVD.



