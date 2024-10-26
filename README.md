[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.07.28-b31b1b.svg?style=flat-square)](https://www.biorxiv.org/content/10.1101/2023.09.05.556135v2)

> **IMPORTANT**  
> This repo is under construction; Docstrings are still being added. Please check back soon!

# A Single-Cell Atlas of ABCA7 Loss-of-Function

This repository contains the main analysis code and links to raw and processed datasets, and relevant analysis pipelines, to reproduce (or extend on) results from our paper.

## Quickstart
- ### Clone the repository and install required packages:
    ```bash
    git clone git@github.com:djunamay/ABCA7lof2.git
    cd ABCA7lof2
    pip install -r requirements.txt
    ```
    - Ensure `mksquashfs` version 4.5 is installed on your system.

- ### See Data Availability below for how to access the data.

- ### See Code Availability for a guide on reproducing the results.

## Data Availability

| Origin | Data Type                          | Raw          | Processed    | Figure Data            | Interactive                                      | Readme |
|--------|------------------------------------|--------------|--------------|------------------------|--------------------------------------------------|--------|
| 🟦 Human | snRNAseq                           | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [dryad](add/link)         | [Broad Single Cell Portal](https://cells.ucsc.edu/)          | [Click](#) |
| 🟦 Human | metadata                           | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [dryad](add/link)         | N/A                                              | [Click](#) |
| 🟦 Human | whole genome sequencing | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [dryad](add/link) | N/A | [wgs_readme.md](data_readmes/wgs_readme.md) |
| 🟦 Human | LC-MS                              | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [syn53461705](https://www.synapse.org/#!Synapse:syn53461705) | [dryad](add/link)         | [MetaboLights](https://www.ebi.ac.uk/metabolights/index)          | [Click](#) |
| 🟩  iPSC | LC-MS                              | [dryad](add/link) | [dryad](add/link) | [dryad](add/link)         | [MetaboLights](https://www.ebi.ac.uk/metabolights/index)           | [Click](#) |
| 🟩  iPSC | confocal images                    | [dryad](add/link) | [dryad](add/link) | [dryad](add/link)         | N/A                                              | [Click](#) |
| 🟩  iPSC | oxygen consumption rates           | [dryad](add/link) | [dryad](add/link) | [dryad](add/link)         | N/A                                              | [Click](#) |
| 🟩  iPSC | biochemical assays                 | [dryad](add/link) | [dryad](add/link) | [dryad](add/link)         | N/A                                              | [Click](#) |

> **N.B.**  
> Regarding data deposited on Synapse ([syn53461705](https://www.synapse.org/#!Synapse:syn53461705)). These data are subject to controlled access in compliance with human privacy regulations. To obtain the data, a data use agreement (DUA) must be completed. This requirement ensures the anonymity of ROSMAP study participants. A DUA can be established with either the Rush University Medical Center (RUMC) or SAGE, the organization that manages Synapse. The necessary forms are available for download on their respective websites.

## Code Availability

| Data Type                | Analysis Type                    | Repository Link                                           |
|--------------------------|-----------------------------------|-----------------------------------------------------------|
| General Analysis         | Replicate Analyses and Figures    | [Run All Analyses](#run-all-analyses)                     |
| Whole Genome Sequencing  | Access Whole-Genome Sequencing    | [ROSMAPwgs](https://github.com/djunamay/ROSMAPwgs)        |
| Genetic Variant Analysis | Get single cell BAM files         | [Get single cell BAM files](#get-single-cell-bam-files)   |
| ABCA7 snRNAseq Analysis  | Perform Celltype Annotation & QC  | [ABCA7 loss-of-function snRNAseq analysis](#abca7-loss-of-function-snrnaseq-analysis) |
| Gene and Pathway Stats   | Compute gene scores and pathway enrichments | [Gene and Pathway Statistics](#gene-and-pathway-statistics) |
| Gene-Pathway             | Perform Gene-Pathway Clustering   | [geneclusters](https://github.com/djunamay/geneclusters)  |
| Confocal Imaging         | Process Confocal Images           | [confocalQuant](https://github.com/djunamay/confocalQuant) |
## Genetic variant analysis

- #### Get single cell BAM files:
    > **<u>Data</u>**  
    > Download the FASTQ files from [Synapse](https://www.synapse.org/#!Synapse:syn53461705).
    - Make the squash file system:
        ```bash
        mksquashfs */fastqs/10x-4819F batch_4819F.sqsh
        mksquashfs */fastqs/10x-4826F batch_4826F.sqsh
        mksquashfs */fastqs/171013Tsa 171013Tsa.sqsh
        ```
    - Run cellranger counting:
        ```bash
        sbatch --array 1-42 */bash_files/cellranger_count.sh
        */bash_files/check_success.sh # iterate over all logs and check whether the pipeline was successful before moving to aggregation
        ```
- #### Run sample swap analysis:
    > **<u>Data</u>**  
    > Access the VCF files..
    - Remap the VCF file:
        ```bash
        /bash_files/crossmap.sh
        */htslib-1.10.2/bgzip out.hg38.vcf --threads 20
        */bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz
        */htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz
        */bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
        */htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz
        ```
    - Run `sample_swap_make_exec.ipynb` to make the text file to iterate through.
    - Run `sample_swap.sh`.
    - Run `./02-sample_swap.ipynb` to visualize sample swap results.

## ABCA7 loss-of-function snRNAseq analysis
> **<u>Data</u>**  
> - Follow [Get single cell BAM files](#get-single-cell-bam-files) above.
> - Then aggregate all the count files by running `./03-aggregate.ipynb`.
> - Alternatively, access the raw aggregated counts matrix, rowData, and colData [here](https://www.synapse.org/#!Synapse:syn53461705).

- #### To Perform Celltype Annotation & QC:
    - Run `./04-get_marker_genes.ipynb` to get marker genes.
    - Run `./05-single_cell_qc_anno.ipynb` to run celltype quality control and annotation.
    - Run `./06-umaps.ipynb` to generate UMAPs.
    - Run `./07-make_sce.ipynb` to save single cell data as SingleCellExperiment object.

- #### To perform ABCA7 loss-of-function-associated statistical analyses:
    - Run `./09-stats_inputs.ipynb` to format data for input to stats analysis.

## Gene and Pathway Statistics
> **<u>Data</u>**  
> - Follow [snRNAseq analysis](#snrna-analysis) above.
> - Alternatively, access the QC'ed aggregated counts matrix, rowData, and colData [here](https://www.synapse.org/#!Synapse:syn53461705). 
> - Or, access the stats_input_data_0825.rds [here](https://figshare.com/s/c944697d9ec30ab06621).

- #### Compute gene scores and pathway enrichments:
    - Run `./10-compute_stats.ipynb`.

- #### Gene score dimensionality reduction and clustering:
    - Run `./11-projections.ipynb`.

- #### K/L Gene Clustering
    - Run `./08-run_partitioning.py` to run METIS and K/L algorithms.
    - Run `./08-benchmarking_graph_partitioning.ipynb` to benchmark clustering and partitioning methods.
    - Run `./12-KL_clusters.ipynb` to visualize graph partitioning results.

## ABCA7 Ala1527Gly analysis
- #### snRNAseq analyses
    > **<u>Data</u>**  
    > - Follow 
    - Run `./19-common_variant_analysis.ipynb`.
    - Run `./20-common_var_plotting.ipynb`.
- #### MD simulations
    > **<u>Data</u>**  

## LC-MS Analyses
- #### Human post-mortem
    > **<u>Data</u>**  
    - Run `./13-plotting_inputs.ipynb`.
    - Run `./16-lipidomics_by_subclass.ipynb` to plot lipidomics aggregate data for the postmortem brain.
- #### iPSC
    > **<u>Data</u>**  
    > - Follow [LC-MS](#lc-ms-analyses) above.
    - Run `./13-plotting_inputs.ipynb`. to get the lipidomic input object.
    - Run `./24-lipidomics_by_subclass.ipynb` to plot lipidomics aggregate data for the iPSC-neuron.
    - Run `./14-figures.ipynb` to plot figures.
    - Run `./25-metabolomics-iN.ipynb` to analyze iPSC-neuron metabolomics data.

## Oxygen Consumption Rates
    > **<u>Data</u>**  
    - Run `./23-seahorse.ipynb`.
    
## Confocal Imaging
> **<u>Data</u>**  
> - Visit the [confocalQuant GitHub Repository](https://github.com/djunamay/confocalQuant).

## Plotting Figures
> **<u>Data</u>**  
    - Run `./14-figures.ipynb` to plot main figure panels.
    - Run `./15-extended-figures.ipynb` to plot extended figures.    
    - Run `./17-variant_carrier_pie_charts.ipynb` to plot variant carrier proportions.
    - Run `./21-basic_pie_charts.ipynb` to plot the common variant analysis.

## Citation

Please cite our preprint using the following BibTeX entry if you use this code in your work:

```bibtex
@article{vonMaydell2023,
  doi = {10.1101/2023.09.05.556135},
  url = {https://doi.org/10.1101/2023.09.05.556135},
  year = {2023},
  month = sep,
  publisher = {Cold Spring Harbor Laboratory},
  author = {Djuna von Maydell and Shannon Wright and Julia Maeve Bonner and Ping-Chieh Pao and Gloria Suella Menchaca and Gwyneth Welch and Carles A. Boix and Hansruedi Mathys and Guillaume Leclerc and Noelle Leary and George Samaan and Manolis Kellis and Li-Huei Tsai},
  title = {A single-cell atlas of {ABCA}7 loss-of-function reveals lipid disruptions, mitochondrial dysfunction and {DNA} damage in neurons}
}
```
