
# ABCA7 variants impact phosphatidylcholine and mitochondria in neurons
`DOI:` [10.1038/s41586-025-09520-y](https://www.nature.com/articles/s41586-025-09520-y)

This repository contains the main analysis code and links to raw and processed datasets, and relevant analysis pipelines, to reproduce (or extend on) results from our paper.

## Setup
- Clone the repository:
```bash
git clone git@github.com:djunamay/ABCA7lof2.git
cd ABCA7lof2
```
  
## Data and Code Availability
- Raw and processed data are available via [Synapse](https://www.synapse.org/#!Synapse:syn53461705), [OSF](https://osf.io/pqr9m/), or [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299277).
- The single-cell atlas can be interactively explored through the [Broad Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP3182).

- See this [table](https://github.com/djunamay/ABCA7lof2/blob/main/open_data/README.md) for a specific overview of data availability. 
- Follow the READMEs below to get to specific code.

| Analysis                                          | Description                                                                                                                                                                     | Go to Code and Readme                                           |
|-----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| Human postmortem - snRNAseq processing          | Process raw fastq files from single-nucleus RNA sequencing of human postmortem samples, including quality control and preprocessing steps.                   | [→](analyses/snRNAseq_processing/)                       |
| Human postmortem - sample swap analysis               | Perform sample-swap checks.                                         | [→](analyses/sample_swap/)                            |
| Human postmortem - snRNAseq stats               | Compute gene- and pathway-level statistics from snRNAseq data.                                         | [→](analyses/snRNAseq_stats/)                            |
| Human postmortem - snRNAseq gene-pathway partitioning  | Evaluate and apply graph partitioning algorithms to gene and pathway-level statistics from snRNAseq data. | [→](analyses/snRNAseq_score_partitioning/)               |
| Human postmortem - snRNAseq ala1527gly          | Analysze snRNAseq data from ROSMAP study participants with the ALA1527GLY variant.          | [→](analyses/common_variant_analysis/)                        |
| induced-neurons - LCMS                            | Analyze lipidomic profiles from iPSC-derived neurons using Liquid Chromatography-Mass Spectrometry (LCMS).             | [→](analyses/iN_LCMS/)                                   |
| induced-neurons - mRNA             |  Analyze mRNA datasets from iPSC-derived neurons.                                  | [→](analyses/bulkRNAseq/) |
| induced-neurons - O2 consumption rates             | Analyze oxygen consumption rates in iPSC neurons            | [→](analyses/iN_O2_consumption/)                         |
| induced-neurons / cortical organoid - imaging              | Analyze MitoHealth, TMRM, and CellROX dyes, and visualize neuronal markers.          | [→](analyses/iN_membrane_potential/)                     |
| cortical organoid - ephys and amyloid               | Plot Amyloid ELISAs and Ephysiological data from cortical organoids.          | [→](analyses/iN_membrane_potential/)                     |

## Software & Package Versions

To install specific environments:

```bash
conda create --name myenv --file env.txt
```
- `bulkrna_packages.txt` and `trim_env_packages.txt` and `qtltools_env_packages.txt` $\rightarrow$ conda environment for bulk mRNA-seq analysis
- `wgs_env_packages.txt` $\rightarrow$ conda environment to access WGS data 
- `confocalquant_packages.txt` $\rightarrow$ conda environment for confocal image analysis
- `bulkrna_packages.txt` $\rightarrow$ conda environment for bulk mRNA-seq analysis
- `scmod_py_packages.txt` $\rightarrow$ conda environment for snRNAseq analysis
- More Info [→](./package_info.md)
