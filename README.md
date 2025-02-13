[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.07.28-b31b1b.svg?style=flat-square)](https://www.biorxiv.org/content/10.1101/2023.09.05.556135v2)

> **IMPORTANT**  
> This repo is under construction; Docstrings are still being added. Please check back soon!

# A Single-Cell Atlas of ABCA7 Loss-of-Function

This repository contains the main analysis code and links to raw and processed datasets, and relevant analysis pipelines, to reproduce (or extend on) results from our paper.

## Setup
- Clone the repository and install required packages:
```bash
git clone git@github.com:djunamay/ABCA7lof2.git
cd ABCA7lof2
pip install -r requirements.txt
```
- Ensure `mksquashfs` version 4.5 is installed on your system.

## Data and Code Availability
| Data                                          | Description                                                                                                                                                                     | Go to Code and Readme                                           |
|-----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| Human postmortem - snRNAseq processing          | Process raw fastq files from single-nucleus RNA sequencing of human postmortem samples, including quality control and preprocessing steps.                   | [$\rightarrow$](analyses/snRNAseq_processing/README.md)                       |
| Human postmortem - snRNAseq stats               | Compute gene- and pathway-level statistics from snRNAseq data.                                         | [$\rightarrow$](analyses/snRNAseq_stats/README.md)                            |
| Human postmortem - snRNAseq score partitioning  | Evaluate and apply graph partitioning algorithms to gene and pathway-level statistics from snRNAseq data. | [$\rightarrow$](analyses/snRNAseq_score_partitioning/README.md)               |
| Human postmortem - snRNAseq ala1527gly          | Analysze snnRNAseq data from ROSMAP study participants with the ALA1527GLY variant.          | [$\rightarrow$](analyses/snRNAseq_ala1527gly/README.md)                        |
| NA - molecular dynamics simulations           | Perform molecular dynamics simulations of ABCA7 protein. | [$\rightarrow$](analyses/molecular_dynamics_simulations/README.md)            |
| iPSC neuron - LCMS                            | Analyze lipidomic profiles from iPSC-derived neurons using Liquid Chromatography-Mass Spectrometry (LCMS).             | [$\rightarrow$](analyses/iN_LCMS/README.md)                                   |
| iPSC neuron - O2 consumption rates             | Analyze oxygen consumption rates in iPSC neurons            | [$\rightarrow$](analyses/iN_O2_consumption/README.md)                         |
| iPSC neuron - membrane potential               | Analyze membrane potential in iPSC neurons via confocal microscopy of membrane potential - dependent dyes           | [$\rightarrow$](analyses/iN_membrane_potential/README.md)                     |
| neurospheroid - calcium dynamics                | Analyze calcium dynamics in neurospheroids to track changes in intracellular calcium levels over time.         | [$\rightarrow$](analyses/neurospheroid_calcium_dynamics/README.md)             |

## Citation

Please cite our preprint using the following BibTeX entry if you use this code in your work:
