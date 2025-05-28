
# snRNAseq Statistics
Compute DEGs and other summary statistics needed for downstream analyses.

## Data availability
- These notebooks take as input the SingleCellExperiment object output by [`analyses/snRNAseq_processing/make_sce.ipynb`](https://github.com/djunamay/ABCA7lof2/blob/main/analyses/snRNAseq_processing/make_sce.ipynb)
  
## Code overview
- Run `./stats_inputs.ipynb` to format data for input to stats analysis.
- Run `./compute_stats.ipynb` to compute statistics used for downstream analyses.
