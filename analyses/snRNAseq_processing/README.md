
## Experiment Description

## Data availability

- FastQ files are deposited on Synapse ([syn53461705](https://www.synapse.org/#!Synapse:syn53461705)). These data are subject to controlled access in compliance with human privacy regulations. To obtain the data, a data use agreement (DUA) must be completed. This requirement ensures the anonymity of ROSMAP study participants. A DUA can be established with either the Rush University Medical Center (RUMC) or SAGE, the organization that manages Synapse. The necessary forms are available for download on their respective websites.
- raw aggregated counts matrix, rowData, and colData [syn53461705](https://www.synapse.org/#!Synapse:syn53461705).
- QC'ed aggregated counts matrix, rowData, and colData [syn53461705](https://www.synapse.org/#!Synapse:syn53461705). 
- VCF files are accessible through Synapse (see this [Git repo](https://github.com/djunamay/ROSMAPwgs) for instructions)
- Where needed, commands to download data are provided at the start of each notebook.

## Code overview

### Run cellranger counting:
- Download the FASTQ files
-  Make the squash file system:
```bash
mksquashfs */fastqs/10x-4819F batch_4819F.sqsh
mksquashfs */fastqs/10x-4826F batch_4826F.sqsh
mksquashfs */fastqs/171013Tsa 171013Tsa.sqsh
```
- Run `../bash_files/cellranger_count.sh`
- Run `../bash_files/check_success.sh` # iterate over all logs and check whether the pipeline was successful before moving to aggregation

### Get metadata files and aggregate counts:
- Run `./metadata.ipynb`
- Run `./aggregate.ipynb`

### Run sample swap analysis:
- Download the VCF files (see this Git repo for instructions)
- Remap the VCF file:
```bash
/bash_files/crossmap.sh
*/htslib-1.10.2/bgzip out.hg38.vcf --threads 20
*/bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz
*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz
*/bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz
```
- Run `../bash_files/sample_swap_make_exec.ipynb` to make the text file to iterate through.
- Run `../bash_files/sample_swap.sh`.
- Run `../bash_files/sample_swap.ipynb` 

### ABCA7 loss-of-function snRNAseq QC & Annotation
- Run `./get_marker_genes.ipynb` to get marker genes.
- Run `./single_cell_qc_anno.ipynb` to run celltype quality control and annotation.
- Run `./umaps.ipynb` to generate UMAPs.
- Run `./make_sce.ipynb` to save single cell data as SingleCellExperiment object.