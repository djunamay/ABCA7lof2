[![bioRxiv](https://img.shields.io/badge/bioRxiv-202023.09.05-b31b1b.svg?style=flat-square)](https://www.biorxiv.org/content/10.1101/2023.09.05.556135v1)

> [!IMPORTANT]  
> This repo is under construction; Docstrings are still being added. Please check back soon!

# A single-cell atlas of ABCA7 loss-of-function 

This repository contains the main analysis code and links to raw and processed datasets, and relevant analysis pipelines, to reproduce (or extend on) results from our paper.  

## Data Availability

All postmortem human data can be accessed through the Synapse AD Knowledge Portal ([syn53461705]((https://www.synapse.org/#!Synapse:syn53461705))), which also includes associated ROSMAP metadata. These data are subject to controlled access in compliance with human privacy regulations. To obtain the data, a data use agreement (DUA) must be completed. This requirement ensures the anonymity of ROSMAP study participants. A DUA can be established with either the Rush University Medical Center (RUMC) or SAGE, the organization that manages Synapse. The necessary forms are available for download on their respective websites. All iPSC-related data are accessible through links provided in our code repositories. For a complete list of data availability and download links, please refer to the code repositories listed below. Additionally, relevant processed datasets are available in the supplementary files of this manuscript.

Follow these instructions to access the data generated and used as part of this study.

- For the processed **WGS data**, follow instructions in the [ROSMAPwgs](https://github.com/djunamay/ROSMAPwgs) repository

> <details>
> <summary>Then run the main.py script with the following parameters</summary>
>
> ```bash
> python main.py --outdir './raw_data/ROSMAP_WGS' --username <USERNAME> --pw <PASSWORD> --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extension 'recalibrated_variants.vcf.gz' --extract_HIGHandMED_annotations False --download True
> python main.py --outdir './raw_data/ROSMAP_WGS' --username <USERNAME> --pw <PASSWORD> --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extension 'annotated.coding.txt' --extract_HIGHandMED_annotations False --download True
> python main.py --outdir './raw_data/ROSMAP_WGS' --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extract_HIGHandMED_annotations True --download False
> ```
> </details>
>
- For the raw and processed **post-mortem snRNAseq data**, **post-mortem lipidomic data**, and patient **metadata** go to [Synapse](https://www.synapse.org/#!Synapse:syn53461705) to request these controlled-access data

- For the raw and processed and raw **neuronal iPSC imaging data**, go to the "Example Use-Case" in [this Github Repository](https://github.com/djunamay/confocalQuant?tab=readme-ov-file#example-use-case) for instructions to download

- For other data files necessary to recapitulate analyses see the respective "Input Data" tabs in the "Analysis" section

For curated data uploades also see:

- All Lipidomic and metabolomic datasets will also be available through the [MetaboLights](https://www.ebi.ac.uk/metabolights/index) database [^1]

- snRNAseq data can be explored on the [USCS Single Cell Browser] soon [^1]


[^1]: Please note that some aspects of these data will be retracted to comply with controlled access regulations. For full access to these data, please visit the repository on Synapse. 

## Code Availability

All code used in this study is available on GitHub. This includes code to replicate the analyses and figure panels presented in the paper [[here]](##analyses), descriptions and code for accessing whole-genome sequencing data [[here]](https://github.com/djunamay/ROSMAPwgs), code for performing gene-pathway clustering [[here]](https://github.com/djunamay/geneclusters), and code for processing confocal images [[here]](https://github.com/djunamay/confocalQuant).

## Methods

For detailed methods descriptions (experimental and computational) please see our paper.

## Run all analyses

#### <u>quickstart</u>

- clone the repository and install required packages

```bash
git clone git@github.com:djunamay/ABCA7lof2.git
pip install -r requirements.txt
```

- follow the steps below as needed

#### <u>single-cell-related</u>

#### To run cellranger counting:
> <details>
> <summary>Input Data</summary>
>
> [Download FASTQ files here](https://www.synapse.org/#!Synapse:syn53461705)    
> </details>
>
> <details>
> <summary>1. make the squash file system</summary>
>
> ```bash
> # Make the squash file systems 
> mksquashfs */fastqs/10x-4819F batch_4819F.sqsh # or modify the cellranger_count.sh script to run without the squash file system
> mksquashfs */fastqs/10x-4826F batch_4826F.sqsh
> mksquashfs */fastqs/171013Tsa 171013Tsa.sqsh
> ```
> </details>
>
> <details>
> <summary>2. run cellranger counting</summary>
>
> ```bash
> # count the FASTQ files:
> sbatch --array 1-42 */bash_files/cellranger_count.sh
> */bash_files/check_success.sh # iterate over all logs and check whether pipeline was successful before moving to aggregation
> ```
> </details>
>
> 3. run *`./03-aggregate.ipynb`* to aggregate all the count files

#### To run sample swap analysis:
>
> <details>
> <summary>Input Data</summary>
>
> See sections **`To run cellranger counting:`**  and **`Data Availability`** above to get BAM files and WGS data.
> </details>
>
> 1. run *`/bash_files/crossmap.sh`* to remap the vcf file
>
> <details>
> <summary>2.then run the following commands in bash:</summary>
> 
> ```bash
> */htslib-1.10.2/bgzip out.hg38.vcf --threads 20 # compress with bgzip
> */bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz # sort the vcf file 
> */htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz # then generate the corresponding tabix file 
> */bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
> *htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz # then generate the corresponding tabix file 
>```
> </details>
>
> 3. run *`sample_swap_make_exec.ipynb`* to make text file to iterate through 
> 4. run *`sample_swap.sh`* to run sample swap 
> 5. run *`./02-sample_swap.ipynb`* to visualize sample swap results 

#### To perform celltype annotation & QC:
> <details>
> <summary>Input Data</summary>
>
> [Download the full aggregated counts matrix, rowData, and colData here](https://www.synapse.org/#!Synapse:syn53461705)    
> </details>
>
> 1. run *`./04-get_marker_genes.ipynb`* to get marker genes for celltype annotation 
> 2. run *`./05-single_cell_qc_anno.ipynb`* to run celltype quality control and annotation I 
> 3. run *`./06-umaps.ipynb`* to run celltype quality control and annotation II 
> 4. run *`./07-make_sce.ipynb`* to save single cell data as singlecellexperiment object 
    
#### To get K/L gene clusters:
> <details>
> <summary>Input Data</summary>
>
> [Download all the necessary input data through figshare] (coming soon)
> </details>
>
> 1. run *`./08-run_partitioning.py`* to run METIS and K\L algorithms 
> 2. run *`./08-benchmarking_graph_partitioning.ipynb`* to benchmark clustering and partitioning methods

#### To perform statistical analyses:
> <details>
> <summary>Input Data</summary>
>
> [Download the full annotated and QCed counts matrix, rowData, and colData here](https://www.synapse.org/#!Synapse:syn53461705) 
>
> [Or download redacted versions (censored patient metadata) through the UCSC Single Cell Browser (coming soon)]
> </details>
>
> 1. run *`./09-stats_inputs.ipynb`* to format data for input to stats analysis 
> 2. run *`./10-compute_stats.ipynb`* to compute gene scores and pathway enrichments 
> 3. run *`./11-projections.ipynb`* for gene score dimensionality reduction and clustering 
> 4. run *`./13-plotting_inputs.ipynb`* to format some data for plotting 
> 5. run *`./19-common_variant_analysis.ipynb`* to compute DEGs for the common ABCA7 variant

#### <u>ipsc-neuron-related</u>

#### For all iPSC Neuronal Omics Analyses:
> <details>
> <summary>Input Data</summary>
> [Download all the necessary input data through figshare] (coming soon)
> </details>  
>
> - see *`./23-seahorse.ipynb`* to visualize graph partitioning results 
> - see *`./24-lipidomics-iN.ipynb`* to plot main figure panels 
> - see *`./25-metabolomics-iN.ipynb`* to plot extended figures

#### For iPSC Neuronal Image Analyses:

- Please go to this [Github Repository](https://github.com/djunamay/confocalQuant)
    
#### <u>plotting/visualization notebooks</u>
> <details>
> <summary>Input Data</summary>
> [Download all the necessary input data through figshare] (coming soon)
> </details>  
>
> - see *`./12-KL_clusters.ipynb`* to visualize graph partitioning results 
> - see *`./14-figures.ipynb`* to plot main figure panels 
> - see *`./15-extended-figures.ipynb`* to plot extended figures
> - see *`./16-lipidomics_PM.ipynb`* to plot lipidomics aggregate data
> - see *`./17-variant_carrier_pie_charts.ipynb`* to plot variant carrier proportions
> - see *`./18-specific_pathway_analysis.ipynb`* to plot genes and pathways for targeted pathway analysis
> - see *`./20-common_var_plotting.ipynb`* to plot the common variant analysis
> - see *`./21-basic_pie_charts.ipynb`* to plot the common variant analysis
> - see *`./22-beta_ox_genes.ipynb`* to plot the common variant analysis


## Citation
Please cite our preprint using the following BibTeX entry if you use this code in your work:
```
@article{vonMaydell2023,
  doi = {10.1101/2023.09.05.556135},
  url = {https://doi.org/10.1101/2023.09.05.556135},
  year = {2023},
  month = sep,
  publisher = {Cold Spring Harbor Laboratory},
  author = {Djuna von Maydell and Shannon Wright and Julia Maeve Bonner and Ping-Chieh Pao and Gloria Suella Menchaca and Gwyneth Welch and Carles A. Boix and Hansruedi Mathys and Guillaume Leclerc and Noelle Leary and George Samaan and Manolis Kellis and Li-Huei Tsai},
  title = {A single-cell atlas of {ABCA}7 loss-of-function reveals lipid disruptions,  mitochondrial dysfunction and {DNA} damage in neurons}
}
```

    
