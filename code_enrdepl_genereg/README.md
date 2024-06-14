## General 
This repository contains scripts for performing enrichment analysis using eFORGE as well as the manual construction of a background CpG set of EPIC array CpGs using the sliding-window approach to match CpGs by methylation levels and variability.

## eFORGE analysis
eFORGE analysis

To run the eFORGE enrichment/depletion analysis, set the working directory to the same directory as where the eFORGE software resides. Next, 
execute the .pl script as follows (example: analysing enrichment in Roadmap DHS hypersensitivity sites)
```
perl eforge.pl -f ./sentinels_194_ID.csv  --array 850k -data erc2-DHS -label erc2-DHS --noproxy 
```

* `2_eforge.pl` : main script for eFORGE enrichment/depletion analysis 
* `2_plot_eforge_enrdepl.R`: calculates permutation p value for assessing enrichment/depletion and generates a heatmap (log2 fold change) for 
