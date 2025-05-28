# snRNA-seq cluster analyses
Partition gene-pathway graph using K/L heuristic.

## Data availability
- Where needed, commands to download data are provided at the start of each notebook.
- This subsection is based on outputs from [`analyses/snRNAseq_stats/compute_stats.ipynb`](https://github.com/djunamay/ABCA7lof2/blob/main/analyses/snRNAseq_stats/compute_stats.ipynb).

## Code overview
- Run `./projections.ipynb` to project DEGs and test for enrichment.
- Run `./run_partitioning.py` to run METIS and K/L algorithms.
- Run `./benchmarking_graph_partitioning.ipynb` to benchmark clustering and partitioning methods.
- Run `./KL_clusters.ipynb` to visualize graph partitioning results.
