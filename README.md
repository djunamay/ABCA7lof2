
# ABCA7 Loss-of-Function Variants Impact Phosphatidylcholine Metabolism in the Human Brain

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
| Human postmortem - snRNAseq processing          | Process raw fastq files from single-nucleus RNA sequencing of human postmortem samples, including quality control and preprocessing steps.                   | [→](analyses/snRNAseq_processing/)                       |
| Human postmortem - snRNAseq stats               | Compute gene- and pathway-level statistics from snRNAseq data.                                         | [→](analyses/snRNAseq_stats/)                            |
| Human postmortem - snRNAseq score partitioning  | Evaluate and apply graph partitioning algorithms to gene and pathway-level statistics from snRNAseq data. | [→](analyses/snRNAseq_score_partitioning/)               |
| Human postmortem - snRNAseq ala1527gly          | Analysze snnRNAseq data from ROSMAP study participants with the ALA1527GLY variant.          | [→](analyses/snRNAseq_ala1527gly/)                        |
| NA - molecular dynamics simulations           | Perform molecular dynamics simulations of ABCA7 protein. | [→](analyses/molecular_dynamics_simulations/)            |
| iPSC neuron - LCMS                            | Analyze lipidomic profiles from iPSC-derived neurons using Liquid Chromatography-Mass Spectrometry (LCMS).             | [→](analyses/iN_LCMS/)                                   |
| iPSC neuron - O2 consumption rates             | Analyze oxygen consumption rates in iPSC neurons            | [→](analyses/iN_O2_consumption/)                         |
| iPSC neuron - membrane potential               | Analyze membrane potential in iPSC neurons via confocal microscopy of membrane potential - dependent dyes           | [→](analyses/iN_membrane_potential/)                     |
| neurospheroid - calcium dynamics                | Analyze calcium dynamics in neurospheroids to track changes in intracellular calcium levels over time.         | [→](analyses/neurospheroid_calcium_dynamics/)             |

