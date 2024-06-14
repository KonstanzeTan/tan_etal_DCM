## Description of individual scripts
eFORGE analysis

To run the eFORGE enrichment/depletion analysis, set the working directory to the same directory as where the eFORGE software resides. Next, 
execute the .pl script as follows (example: analysing enrichment in Roadmap DHS hypersensitivity sites)
```
perl eforge.pl -f ./sentinels_194_ID.csv  --array 850k -data erc2-DHS -label erc2-DHS --noproxy 
```

* `maf` : the minor allele frequency threshold to be applied when converting impute2 genotypes to dosages.
* `post_maf` : the minor allele frequency threshold to be applied after matrixEQTL has been run. useful to define different summaries.
* `cisDist` : the maximal distance between SNP-CpGs pairs until which it counts as a 'cis-pairs'. cis-pairs will be reported in a separate file.
* `pvThresholdCis` : The p-value threshold to be applied for cis-meQTLs. Only associations with `p-value < pvThresholdCis` will be reported
* `pvThresholdTrans` : The p-value threshold to be applied for trans-meQTLs. Only associations with `p-value < pvThresholdTrans` will be reported
* `significant_pv_cutoff` : For final results files: at which threshold should an association be called significant? e.g.: 1e-14
