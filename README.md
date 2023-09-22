***This repository contains code to reproduce the analyses in "ABCA7 loss of function induces lipid dyshomeostasis and DNA damage in neurons"***

### 1. Data Availability

*Follow these instructions to access the data generated and used as part of this study.*

- For whole genome sequencing data
> Follow instructions in the [ROSMAPwgs](https://github.com/djunamay/ROSMAPwgs) repository

- For snRNA-seq and post-mortem lipidomic data
> Download from **Synapse** [coming soon](https://linktosynapse) \
> This is what the file structure on Synapse looks like: \
> click [here](https://linktosynapse) for descriptions
> ```
> ABCA7lof
> └───fastq
>     └───10x-4819F
>     └───10x-4826F
>     └───171013Tsa
>     └───meta
>     └───md5sums
> └───raw_data
>     └───raw_gene_names.csv 
>     └───raw_sample_metadata.csv 
>     └───raw_counts.mtx 
> └───qc_data
>     └───qc_gene_names.csv
>     └───qc_sample_metadata.csv 
>     └───qc_counts.mtx 
> └───pm_lipidomic_data
>     └───
> ```

- For all other data
> Download from the **Open Science Framework**, [here](https://osf.io/vn7w2/) \
> This is what the file structure on OSF looks like: \
> click [here](https://osf.io/vn7w2/) for descriptions
> ```
> ABCA7lof
> └───single_cell_data
>     └───
> └───ipsc_data
>     └───
> └───other_postmortem_data
>     └───
> └───supplementary_tables
>     └───
> ```

### 2. Data Analysis

*If you'd like to start with the fastq files:*

a. **`Download the fastq files` [[here](https://linktosynapse)]** 


b. **`run cellranger counting`** 

> ```bash
> # Make the squash file systems 
> mksquashfs */fastqs/10x-4819F batch_4819F.sqsh # or modify the cellranger_count.sh script to run without the squash file system
> mksquashfs */fastqs/10x-4826F batch_4826F.sqsh
> mksquashfs */fastqs/171013Tsa 171013Tsa.sqsh
> ```

> ```bash
> # count the FASTQ files:
> sbatch --array 1-42 */bash_files/cellranger_count.sh
> */bash_files/check_success.sh # iterate over all logs and check whether pipeline was successful before moving to aggregation
> ```

c. **`sample swap`**

*NB. We did this for the analysis as a control to check that WGS data and snRNA-seq data match (they do), so you don't need to run this again*
> see *`/bash_files/crossmap.sh`* to remap the vcf file \
> then run the following commands in bash:
> ```bash
> */htslib-1.10.2/bgzip out.hg38.vcf --threads 20 # compress with bgzip
> */bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz # sort the vcf file 
> */htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz # then generate the corresponding tabix file 
> */bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
> *htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz # then generate the corresponding tabix file 
>```
> see *`sample_swap_make_exec.ipynb`* to make text file to iterate through \
> see *`sample_swap.sh`* to run sample swap \
> see *`./02-sample_swap.ipynb`* to visualize sample swap results 

d. **`aggregate counts`** 
> see *`./03-aggregate.ipynb`* to aggregate all the count files

e. **`celltype annotation & QC`** 

*If you'd like to start with the annotations, start here:*

- download x here

> see *`./04-get_marker_genes.ipynb`* to get marker genes for celltype annotation \
> see *`./05-single_cell_qc_anno.ipynb`* to run celltype quality control and annotation I \
> see *`./06-umaps.ipynb`* to run celltype quality control and annotation II \
> see *`./07-make_sce.ipynb`* to save single cell data as singlecellexperiment object 

e. **`gene clusters`** 

*If you'd like to play with the gene clustering, start here:*
> run *`./08-run_partitioning.py`* to run METIS and K\L algorithms \
> see *`./08-processing_gsets.ipynb`* to benchmark clustering and partitioning methods

<details>
<summary>Methods</summary>
<br>
methods from paper
</details>


f. **`stats`**

*If you'd like to start with the stats, start here:*

> see *`./09-stats_inputs.ipynb`* to format data for input to stats analysis \
> see *`./10-compute_stats.ipynb`* to compute gene scores and pathway enrichments \ 

g. **`plots`**

*If you'd like to just replot:*
> see *`./11-projections.ipynb`* for gene score dimensionality reduction and clustering
> see *`./12-KL_clusters.ipynb`* to visualize graph partitioning results \
> see *`./13-plotting_inputs.ipynb`* to format some data for plotting \
> see *`./14-figures.ipynb`* to plot main figure panels \
> see *`./15-extended-figures.ipynb`* to plot extended figures

