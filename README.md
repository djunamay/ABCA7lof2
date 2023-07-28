***This repository contains code to reproduce the analyses presented in***
### ABCA7 loss of function induces lipid dyshomeostasis and DNA damage in neurons

#### Data Availability

- Fastq files and metadata
> The raw data repository on Synapse is organized as follows (comments indicate the code that generated these data).
> ```
> ABCA7lof
> └───fastq
>     └───10x-4819F
>     └───10x-4826F
>     └───171013Tsa
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
> mksquashfs */fastqs/10x-4819F batch_4819F.sqsh
> mksquashfs */fastqs/10x-4826F batch_4826F.sqsh
> mksquashfs */fastqs/171013Tsa 171013Tsa.sqsh
> ```

> ```bash
> # count and aggregate the FASTQ files:
> sbatch --array 1-42 */bash_files/cellranger_count.sh
> */bash_files/check_success.sh # iterate over all logs and check whether pipeline was successful before moving to aggregation
> ```

4. Run **`./sample swap`** to perform sample swap analysis

> ```bash
>*/bash_files/crossmap.sh #to remap the vcf file (see crossmap.sh)
>*/htslib-1.10.2/bgzip out.hg38.vcf --threads 20 #compress with bgzip
>*/bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz #sort the vcf file: 
>*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz #then generate the corresponding tabix file 
>*/bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
>*htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz # then generate the corresponding tabix file 
> ```

Follow  **`*/bash_files/sample_swap_make_exec.ipynb`**

> ```bash
> sbatch --array 1-42 */bash_files/sample_swap.sh
> ```

5. Follow **`./aggregate.ipynb`** to aggregate the individual count profiles into single matrix

5. Follow  **`./single_cell_qc_anno.ipynb`** to QC and ANNOTATE the single-cell counts matrix:

6. Follow **`./pathway_analsysis.ipynb`** and **`./other_analysis.ipynb`** to perform all analyses

7. Follow **`./figures.ipynb`** to plot all the figures
