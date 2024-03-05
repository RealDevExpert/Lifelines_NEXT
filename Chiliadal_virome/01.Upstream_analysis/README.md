Directories contain attempts at code optimization and often contain their own README (this feature will probably be removed to avoid confusion).

The scripts outside directories are the actual scripts used for the upstream data analysis and are enumerated in the order of execution (the same numbers are used for scripts that can be executed in parallel). Meaning of some cryptic two-letter abbreviations after the number: 

- sc: single-cell
- vd: virus discovery
- pd: post (virus) discovery
- at: abundance table

To execute table_of_origin.R in the upstream analysis (do that between 06. and 08. steps):

```
# DIR: /scratch/p282752/ANALYSIS_CHILIADAL/scripts
# DIR exists: /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/table_of_origin
module load R
# mind the / at the end of the first argument
Rscript table_of_origin.R /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/table_of_origin/ *_table_of_origin
```

To execute combine_coverage.R in the upstream analysis (do that after the dereplication):

```
module load R
# mind the / at the end of the first argument
Rscript combine_coverage.R /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/mapping/VLP_to_w_neg_der95/coverage/ *.coverage.txt
```
