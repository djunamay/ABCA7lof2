***This repository contains code to reproduce the analyses in "ABCA7 loss of function induces lipid dyshomeostasis and DNA damage in neurons"***

### 1. Data Availability

- WGS data
to get genomic variant data 
> Follow instructions here:
> https://github.com/djunamay/ROSMAPwgs

- Single cell data: download from **Synapse** [here](link to synapse).
> ```
> ABCA7lof
> └───fastq
>     └───10x-4819F
>     └───10x-4826F
>     └───171013Tsa
>     └───meta
>     └───md5sums
> └───raw_data
>     └───raw_gene_names.csv # comment on origin
>     └───raw_sample_metadata.csv 
>     └───raw_counts.mtx 
> └───qc_data
>     └───qc_gene_names.csv
>     └───qc_sample_metadata.csv 
>     └───qc_counts.mtx 
> └───pm_lipidomic_data
>     └───
> ```

- Other data generated as part of this study: download from the **Open Science Framework** [here](https://osf.io/vn7w2/). 

- Other data analyzed in this study


### 2. Data pre-processing
**`./1-metadata.ipynb`**
Extract metadata for post-mortem samples and match it to the snRNAseq library IDs. Extract genomic variant info of interest and match that as well.

**`cellranger counting`** 

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

**`sample swap`**
> ```bash
>*/bash_files/crossmap.sh #to remap the vcf file (see crossmap.sh)
>*/htslib-1.10.2/bgzip out.hg38.vcf --threads 20 #compress with bgzip
>*/bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz #sort the vcf file: 
>*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz #then generate the corresponding tabix file 
>*/bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
>*htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz # then generate the corresponding tabix file 
> `

**`./02-sample_swap.ipynb`**
Visualize sample swap results

**`./03-aggregate.ipynb`**
Aggregate cellranger count outputs

**`./04-get_marker_genes.ipynb`**
Get marker genes for celltype annotation

**`./05-single_cell_qc_anno.ipynb`**
Run celltype quality control and annotation I

**`./06-umaps.ipynb`**
Run celltype quality control and annotation II

**`./07-make_sce.ipynb`**
Save single cell data as singlecellexperiment object

### 3. Stats
**`./08-processing_gsets.ipynb`**
Processing genesets and evaluating KL heuristic

**`./09-stats_inputs.ipynb`**
Format data for input to stats analysis

**`./10-compute_stats.ipynb`**
Compute DEGs and pathway enrichments 

**`./11-projections.ipynb`**
Projecting DEGs into UMAP space and clustering

**`./12-KL_clusters.ipynb`**
Assign leading edge genes to clusters based on gene-pathway graph


### 4. Figures
**`./13-plotting_inputs.ipynb`**
Process some of the external/non-single cell datasets for plotting

**`./14-figures.ipynb`**
Plot remaining main figure panels

**`./15-extended-figures.ipynb`**
Plot extended figures

