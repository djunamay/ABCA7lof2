***This repository contains code to reproduce the analyses presented in***
### ABCA7 loss of function induces lipid dyshomeostasis and DNA damage in neurons

#### Data Availability

- Fastq files and metadata
> The raw data repository on Synapse is organized as follows (comments indicate the code that generated these data).
> ```
> ABCA7lof
> └───fastq
>     └───fastq1..
>     └───...
>     └───meta
>     └───md5sums
> ```

- If you would like to *process the raw counts matrix* and associated metadata, these files can be downloaded on **Synapse** [here](link to synapse). [coming soon]
> The raw data repository on Synapse is organized as follows (comments indicate the code that generated these data).
> ```
> ABCA7lof
> └───raw_data
>     └───raw_gene_names.csv # comment on origin
>     └───raw_sample_metadata.csv 
>     └───raw_counts.mtx 
> ```

- If you would like to *access the fully-processed, annotated, and qc-ed data*, that data can be found on **Synapse** [here](link to synapse). [coming soon]
> The qc-ed data repository on Synapse is organized as follows (comments indicate the code that generated these data).
> ```
> ABCA7lof
> └───qc_data
>     └───qc_gene_names.csv
>     └───qc_sample_metadata.csv 
>     └───qc_counts.mtx 
> ```

- All other data required to reproduce cell type qc, annotation, and analyses are available on the **Open Science Framework** [here](https://osf.io/vn7w2/).
> The data repository on OSF is organized as follows (comments indicate the code that generated these data).
> ```

> ```

- External datasets used for analysis:

    **Proteomic Dataset:** https://www.synapse.org/#!Synapse:syn21449447
    
#### 

1. Run **`./wgs script`** to get genomic variant data 
> Follow instructions here:
> https://github.com/djunamay/ROSMAPwgs
>

2. Run **`./metadata.ipynb`** to get all the necessary metadata

3. Run **`./bash_files/cellranger_count.sh`** and **`./bash_files/aggregate.sh`** to align, count, and aggregate the FASTQ files.

> ```bash
> # Make the squash file systems 
> mksquashfs ../../data/fastqs/10x-4819F batch_4819F.sqsh
> mksquashfs ../../data/fastqs/10x-4826F batch_4826F.sqsh
> mksquashfs ../../data/fastqs/171013Tsa 171013Tsa.sqsh

> # count and aggregate the FASTQ files:
> Rscript ./preprocessing_scripts/counting/get_cellranger_count_exec.r
> sbatch ./preprocessing_scripts/counting/cellranger_count.sh
> Rscript ./preprocessing_scripts/counting/get_metadata_for_cellranger_aggr.r
> sbatch ./preprocessing_scripts/counting/aggregate.sh
> ```

4. Run **`./sample swap`** to perform sample swap analysis

5. Follow  **`./single_cell_qc_anno.ipynb`** to QC and ANNOTATE the single-cell counts matrix:

6. Follow **`./pathway_analsysis.ipynb`** and **`./other_analysis.ipynb`** to perform all analyses

7. Follow **`./figures.ipynb`** to plot all the figures
