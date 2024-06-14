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




