
# ABCA7 Loss-of-Function Variants Impact Phosphatidylcholine Metabolism in the Human Brain

This repository contains the main analysis code and links to raw and processed datasets, and relevant analysis pipelines, to reproduce (or extend on) results from our paper.

## Setup
- Clone the repository and install required packages:
```bash
git clone git@github.com:djunamay/ABCA7lof2.git
cd ABCA7lof2
```
- See analysis-specific subfolders below for package requirements, data links, and code
  
## Data and Code Availability
| Analysis                                          | Description                                                                                                                                                                     | Go to Code and Readme                                           |
|-----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| Human postmortem - snRNAseq processing          | Process raw fastq files from single-nucleus RNA sequencing of human postmortem samples, including quality control and preprocessing steps.                   | [→](analyses/snRNAseq_processing/)                       |
| Human postmortem - sample swap analysis               | Perform sample-swap checks.                                         | [→](analyses/sample_swap/)                            |
| Human postmortem - snRNAseq stats               | Compute gene- and pathway-level statistics from snRNAseq data.                                         | [→](analyses/snRNAseq_stats/)                            |
| Human postmortem - snRNAseq gene-pathway partitioning  | Evaluate and apply graph partitioning algorithms to gene and pathway-level statistics from snRNAseq data. | [→](analyses/snRNAseq_score_partitioning/)               |
| Human postmortem - snRNAseq ala1527gly          | Analysze snnRNAseq data from ROSMAP study participants with the ALA1527GLY variant.          | [→](analyses/common_variant_analysis/)                        |
| NA - molecular dynamics simulations           | Perform molecular dynamics simulations of ABCA7 protein. | [→](analyses/molecular_dynamics_simulations/)            |
| induced-neurons - LCMS                            | Analyze lipidomic profiles from iPSC-derived neurons using Liquid Chromatography-Mass Spectrometry (LCMS).             | [→](analyses/iN_LCMS/)                                   |
| induced-neurons - mRNA             |  Analyze mRNA datasets from iPSC-derived neurons.                                  | [→](analyses/iN_LCMS/) |
| induced-neurons - O2 consumption rates             | Analyze oxygen consumption rates in iPSC neurons            | [→](analyses/iN_O2_consumption/)                         |
| induced-neurons - fixed dyes              | Analyze MitoHealth imaging.          | [→](analyses/iN_membrane_potential/)                     |
| induced-neurons - live dyes               | Analyze TMRM and CellROX imaging.           | [→](analyses/iN_membrane_potential/)                     |
| cortical organoid - other plots               | Plot Amyloid ELISAs and Ephysiological data from cortical organoids.          | [→](analyses/iN_membrane_potential/)                     |
| induced-neurons and cortical organoids - marker images               | Visualize neuronal marker imaging.         | [→](analyses/iN_membrane_potential/)                     |
