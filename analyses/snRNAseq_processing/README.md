
# snRNA-seq data processing (QC, annotation)

## Experiment Description

### Isolation of nuclei from frozen postmortem brain tissue

**For batch #1:** The protocol for the isolation of nuclei from frozen postmortem brain tissue (region BA10) was adapted for smaller sample volumes from a previous study [@Mathys2019-wb]. All procedures were carried out on ice or at 4°C. In brief, postmortem brain tissue was homogenized in 700 µl Homogenization Buffer (320 mM sucrose, 5 mM CaCl₂, 3 mM Mg(CH₃COO)₂, 10 mM Tris HCl pH 7.8, 0.1 mM EDTA pH 8.0, 0.1% IGEPAL CA-630, 1 mM β-mercaptoethanol, and 0.4 U/µl recombinant RNase inhibitor (Clontech)) using a Wheaton Dounce tissue grinder (15 strokes with the loose pestle). Homogenized tissue was filtered through a 40-µm cell strainer and mixed with an equal volume of Working Solution, prepared by mixing Diluent (30 mM CaCl₂, 18 mM Mg(CH₃COO)₂, 60 mM Tris pH 7.8, 0.6 mM EDTA, 6 mM β-mercaptoethanol) with Optiprep density gradient solution (Sigma-Aldrich D1556-250ML) in a 1:5 ratio. The sample mix was then layered onto an Optiprep density gradient composed of 750 µl 30% OptiPrep solution (1.5:1 ratio of Working Solution:Homogenization Buffer) over 300 µl 40% OptiPrep solution (4:1 ratio of Working Solution:Homogenization Buffer). The nuclei were separated by centrifugation (5 min, 10,000 × g, 4°C). Approximately 100 µl of nuclei were collected from the 30%/40% interface and washed twice with 1 ml of PBS containing 0.04% BSA, centrifuging at 300 × g for 3 min (4°C) between washes. Nuclei were then resuspended in 100 µl PBS containing 0.04% BSA, counted using a C-Chip disposable hemocytometer, and diluted to 1000 nuclei/µl in PBS containing 0.04% BSA.

**For batch #2:** These samples (fresh postmortem brain; PFC BA10) were prepared as part of and according to a previous study [@Mathys2019-wb].

Informed consent and Anatomical Gift Act consent were obtained from each participant. The Religious Orders Study and Rush Memory and Aging Project were approved by the Institutional Review Board (IRB) of Rush University Medical Center. All participants signed a repository consent, allowing their data and biospecimens to be shared.

### Droplet-based snRNA-seq

**For batch #1:** cDNA libraries were generated using the Chromium Single Cell 3′ Reagent Kits v3 following the manufacturer's protocol (10x Genomics). Libraries were sequenced on the NovaSeq 6000 S2 platform (paired-end, 28 + 91 bp, with an 8-nucleotide index). Samples were distributed across two lanes and sequenced twice on separate flow cells to enhance sequencing depth.

**For batch #2:** Libraries were prepared using Chromium Single Cell 3′ Reagent Kits v2 and sequenced with the NextSeq 500/550 High Output v2 kits (150 cycles), as described in our previously published study [@Mathys2019-wb].

Raw sequencing reads from all samples were processed jointly for alignment and gene counting.

## Data availability

- FastQ files are deposited on Synapse ([syn53461705](https://www.synapse.org/#!Synapse:syn53461705)). These data are subject to controlled access in compliance with human privacy regulations. To obtain the data, a data use agreement (DUA) must be completed. This requirement ensures the anonymity of ROSMAP study participants. A DUA can be established with either the Rush University Medical Center (RUMC) or SAGE, the organization that manages Synapse. The necessary forms are available for download on their respective websites.
- raw aggregated counts matrix, rowData, and colData [syn53461705](https://www.synapse.org/#!Synapse:syn53461705).
- QC'ed aggregated counts matrix, rowData, and colData [syn53461705](https://www.synapse.org/#!Synapse:syn53461705). 
- VCF files are accessible through Synapse (see this [Git repo](https://github.com/djunamay/ROSMAPwgs) for instructions)
- Where needed, commands to download data are provided at the start of each notebook.

## Code overview

### Requirements:
- `mksquashfs` (version 4.5)

### Run cellranger counting:
- Download the FASTQ files (see above)
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
<!-- 
### Run sample swap analysis:
- Download the VCF files (see [this Git repo](https://github.com/djunamay/ROSMAPwgs) for instructions)
- Remap the VCF file:
```bash
../../bash_files/crossmap.sh
*/htslib-1.10.2/bgzip out.hg38.vcf --threads 20
*/bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz
*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz
*/bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz
```
- Run `../../bash_files/sample_swap_make_exec.ipynb` to make the text file to iterate over.
- Run `../../bash_files/sample_swap.sh`.
- Run `./sample_swap.ipynb`  -->

### ABCA7 loss-of-function snRNAseq QC & Annotation
- Run `./get_marker_genes.ipynb` to get marker genes.
- Run `./single_cell_qc_anno.ipynb` to run celltype quality control and annotation.
- Run `./umaps.ipynb` to generate UMAPs.
- Run `./make_sce.ipynb` to save single cell data as SingleCellExperiment object.

