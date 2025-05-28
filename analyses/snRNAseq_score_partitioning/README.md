# snRNA-seq cluster analyses
Project genes into lower-dimensional embedding based on cell type perturbation scores.
Partition gene-pathway graph using K/L heuristic.

## Data availability
- Where needed, commands to download data are provided at the start of each notebook.
- This subsection is based on outputs from [`analyses/snRNAseq_stats/compute_stats.ipynb`](https://github.com/djunamay/ABCA7lof2/blob/main/analyses/snRNAseq_stats/compute_stats.ipynb).

## Code overview
- Run `./projections.ipynb` to project DEGs and test for enrichment.
- Run `./run_partitioning.py` to run METIS and K/L algorithms.
- Run `./benchmarking_graph_partitioning.ipynb` to benchmark clustering and partitioning methods.
- Run `./KL_clusters.ipynb`[^1] to visualize graph partitioning results.

[1]: Plots in the paper were generated using [/analyses/bulkRNAseq/KL_bulk.ipynb](https://github.com/djunamay/ABCA7lof2/blob/main/analyses/bulkRNAseq/KL_bulk.ipynb) to allow for consistent coloring with the iN neuronal clusters.
