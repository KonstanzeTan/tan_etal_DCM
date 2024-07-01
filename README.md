# Code for Epigenome-wide association analysis


## Table of Contents


   * [Quantification and statistical analysis of DNA methylation](#quantification-and-statistical-analysis-of-dna-methylation)
   * [Enrichment analysis of genomic regulatory features](#enrichment-analysis-of-genomic-regulatory-features)
   * [Expression quantitative trait methylation analysis](#expression-quantitative-trait-methylation-analysis)
   * [Causal analyses (Mendelian Randomisation and Colocalisation)](#causal-analyses-mendelian-randomisation-and-colocalisation)
   * [Weighted Gene Correlation Network Analysis (WGCNA)](#weighted-gene-correlation-network-analysis-wgcna)
   * [Construction of Methylation Risk Score (MRS)](#construction-of-methylation-risk-score-mrs)



## Quantification and statistical analysis of DNA methylation 

We ran EWAS of DCM in two independent cohorts to identify conditionally-independent sentinel CpGs (P<5.96E-08) with confirmed consistent directionality of association. Within the multi-ancestry cohort (MAGNet), ancestry-specific EWAS was performed followed by inverse variance-weighted meta-analysis using [METAL](https://csg.sph.umich.edu/abecasis/metal/). A 2-stage correction of inflation in effect sizes and standard errors was applied: once to ancestry-specific results prior to meta-analysis, and once after meta-analysis. Candidate sentinel CpGs were identified (Bonferroni P<0.05), following which sentinel CpGs were identified by assignment into distinct genomic loci as well as conditional analysis. 

The main EWAS model is

```
DCM ~beta(quantile-normalised) + Age + Gender + principal components of control probes capturing ~95% of control probe variation 
```

The scripts in [Quantification and statistical analysis of DNA methylation](./code_quantmethylationstatanalysis/) facilitate various tasks in methylation data preprocessing and quality control, calculation of methylation betas, epigenome-wide association analysis, meta-analysis, adjustment for test-statistic inflation and finally the identification of sentinel CpGs. 

## Enrichment analysis of genomic regulatory features

We analysed sentinel CpGs for enrichment in various genomic regulatory features using a permutation testing approach, whereby sentinel CpG overlap of genomic regulatory features is compared to a background set comprising permutations of other EPIC array CpGs matched to sentinel CpGs by genomic location, methylation levels and/or variability. This addresses bias inherent in methylation arrays, which preferentially assay pre-determined genomic sites and well-annotated genes

Enrichment analysis for tissue-specific chromatin states, histone-marked regions and DNAse 1 hypersensitive sites was conducted using a local installation of eFORGE v2 downloaded from [Altius Institute eForge](https://eforge.altiusinstitute.org/). Besides enrichment, depletion was also investigated by specifying the  `--depletion` flag. Significant enrichment or depletion was not assessed using the default binomial P value output by the software. Instead, we defined significant enrichment or depletion by calculating a permutation P value usign overlap counts from each of the background permutation sets. eFORGE auto-constructs a background set of permuted CpGs matched to sentinel CpGs by CpG island and gene annotation. 

Enrichment analysis for expression quantitative trait methylation loci (eQTM) and transcription factor binding sites (TFBS) were assessed using a manually constructed background set comprising permutations of EPIC array CpGs matched to sentinel CpGs by methylation levels and variability.

The scripts in [Enrichment analysis of genomic regulatory features](./code_enrdepl_genereg/) were used for enrichment analysis and calculation of permutation p values (eFORGE) as well as manual construction of background set of permuted EPIC array CpGs. 

## Expression quantitative trait methylation analysis

We identified significant associations between sentinel CpG methylation and proximal gene expression (+/- 1Mb) of the sentinel CpG. Replicated eQTMs were defined as significant eQTMs in the discovery cohort (MAGNet; disc FDR P<0.05) for which replication testing (BMCB) confirmed consistent directionality of association. Association testing was performed using the R package MatrixEQTL using the following model: 

```
gene expression ~  methylation + Age + Gender + ancestry + RNA Integrity Number(RIN) + 5 PEER factors
```

PEER factors correspond to hidden sources of variation, learned by Bayesian probabilistic estimation of residual factors to accurately capture and account for latent influences in the data, thereby enhancing the robustness of statistical models. The choice of number of PEER factors to learn was informed by the [GTex eQTL mapping study](https//www.nature.com/articles/nature24277): for N < 150, use 15 PEERs; 150 <= N < 250, use 30 PEERs; N >= 250, use 35 PEERs.  

The scripts in [Expression quantitative trait methylation analysis](./code_eqtm/) were used for learning PEER factors from gene expression data and implementing linear models to identify eQTMs (MatrixEQTL).

## Causal analyses (Mendelian Randomisation and Colocalisation) 

We elucidated the putative causal contribution of sentinel CpG methylation to DCM and proximal gene expression using [SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview). Significant associations were validated using one-sample MR and underwent colocalisation analysis to assess the posterior probability of a shared causal variant underlying the association between the two assessed traits 

The primary SMR analysis utilizes SNPs linked to sentinel CpG methylation (known as 'methylation quantitative trait loci') derived from MAGNet left ventricular samples in the current investigation, as well as publicly-available genetic associations from the following sources:
  - [GWAS of DCM in UKBiobank cohort](https://humandbs.biosciencedbc.jp/en/hum0197-v3)
  - [Left ventricular eQTL from GTex v8 release](https://console.cloud.google.com/storage/browser/gtex-resources?pli=1)

Individual-level genotype, methylation and expression/disease status data from the MAGNet cohort were used for one-sample MR analyses.

The scripts in [Causal analyses](./code_causal_analyses/) detail the primary SMR analysis, one-sample MR validation of SMR-significant hits, and colocalisation analysis of genetic associations for the assessed traits. 

## Weighted gene correlation network analysis (WGCNA)

We employed WGCNA to detect co-methylation modules amongst DCM-associated CpGs (discovery EWAS FDR<0.05, directionality of association confirmed in replication cohort, SD>0.02). 

The scripts in [Co-methylation analyses](./code_wgcna/) identify common co-methylation patterns by comparing ancestry-specific co-methylation networks.

## Construction of methylation risk score (MRS)

Following fine-mapping of 28 DCM sentinel CpGs and testing individual CpGs within fine-mapped regions for CVD trait associations, we constructed CVD-trait specific MRS as a weighted sum of methylation betas from CpGs within fine-mapped regions and tested MRS for association with CVD.

The scripts in [MRS](./code_MRS/) detail the construction of region and trait-specific MRS and testing their associations with the respective CVD traits.



