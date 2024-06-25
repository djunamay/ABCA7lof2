[![bioRxiv](https://img.shields.io/badge/bioRxiv-202023.09.05-b31b1b.svg?style=flat-square)](https://www.biorxiv.org/content/10.1101/2023.09.05.556135v1)

> **IMPORTANT**  
> This repo is under construction; Docstrings are still being added. Please check back soon!

# A Single-Cell Atlas of ABCA7 Loss-of-Function

This repository contains the main analysis code and links to raw and processed datasets, and relevant analysis pipelines, to reproduce (or extend on) results from our paper.

## Data Availability

- **Postmortem Human Data**: Accessible through the Synapse AD Knowledge Portal ([syn53461705](https://www.synapse.org/#!Synapse:syn53461705)), which also includes associated ROSMAP metadata[^1].
- **iPSC-Related Data**: Accessible through links provided [below](#ipsc-neuron-related).
- **Processed WGS Data**: Follow instructions in the [ROSMAPwgs](https://github.com/djunamay/ROSMAPwgs) repository.

<details>
<summary>Run the main.py script with the following parameters</summary>

```bash
python main.py --outdir './raw_data/ROSMAP_WGS' --username <USERNAME> --pw <PASSWORD> --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extension 'recalibrated_variants.vcf.gz' --extract_HIGHandMED_annotations False --download True
python main.py --outdir './raw_data/ROSMAP_WGS' --username <USERNAME> --pw <PASSWORD> --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extension 'annotated.coding.txt' --extract_HIGHandMED_annotations False --download True
python main.py --outdir './raw_data/ROSMAP_WGS' --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extract_HIGHandMED_annotations True --download False
```
</details>

- For other data files necessary to recapitulate analyses, see the respective `Input Data` tabs in the [Run All Analyses](#run-all-analyses) section.
- **Lipidomic and Metabolomic Datasets**: Available through the [MetaboLights](https://www.ebi.ac.uk/metabolights/index) database[^2].
- **snRNAseq Data**: Can be explored on the [UCSC Single Cell Browser](https://cells.ucsc.edu/) soon[^2].

[^1]: These data are subject to controlled access in compliance with human privacy regulations. To obtain the data, a data use agreement (DUA) must be completed. This requirement ensures the anonymity of ROSMAP study participants. A DUA can be established with either the Rush University Medical Center (RUMC) or SAGE, the organization that manages Synapse. The necessary forms are available for download on their respective websites.
[^2]: Please note that some aspects of these data will be retracted to comply with controlled access regulations. For full access to these data, please visit the repository on Synapse.

## Code Availability

- [Replicate the analyses and figure panels presented in the paper](#run-all-analyses)
- [Access whole-genome sequencing data](https://github.com/djunamay/ROSMAPwgs)
- [Perform gene-pathway clustering](https://github.com/djunamay/geneclusters)
- [Process confocal images](https://github.com/djunamay/confocalQuant)

## Methods

For detailed methods descriptions (experimental and computational), please see our paper.

## Run All Analyses

### Quickstart

1. **Clone the repository and install required packages**:
    ```bash
    git clone git@github.com:djunamay/ABCA7lof2.git
    cd ABCA7lof2
    pip install -r requirements.txt
    ```

2. **Install other packages**:
    Ensure `mksquashfs` version 4.5 is installed on your system.

3. **Follow the steps below as needed**.

### Single-Cell Related

#### To Run Cellranger Counting:

<details>
<summary>Input Data</summary>
<a href="https://www.synapse.org/#!Synapse:syn53461705">Download FASTQ files here</a>
</details>

1. **Make the squash file system**:
    ```bash
    mksquashfs */fastqs/10x-4819F batch_4819F.sqsh
    mksquashfs */fastqs/10x-4826F batch_4826F.sqsh
    mksquashfs */fastqs/171013Tsa 171013Tsa.sqsh
    ```

2. **Run cellranger counting**:
    ```bash
    sbatch --array 1-42 */bash_files/cellranger_count.sh
    */bash_files/check_success.sh # iterate over all logs and check whether the pipeline was successful before moving to aggregation
    ```

3. **Aggregate all the count files**:
    - Run `./03-aggregate.ipynb`.

#### To Run Sample Swap Analysis:

<details>
<summary>Input Data</summary>
See sections **To Run Cellranger Counting** and **Data Availability** above to get BAM files and WGS data.
</details>

1. **Remap the VCF file**:
    ```bash
    /bash_files/crossmap.sh
    ```

2. **Run the following commands in bash**:
    ```bash
    */htslib-1.10.2/bgzip out.hg38.vcf --threads 20
    */bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz
    */htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz
    */bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
    */htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz
    ```

3. **Run sample swap**:
    - Run `sample_swap_make_exec.ipynb` to make the text file to iterate through.
    - Run `sample_swap.sh`.
    - Run `./02-sample_swap.ipynb` to visualize sample swap results.

#### To Perform Celltype Annotation & QC:

<details>
<summary>Input Data</summary>
<a href="https://www.synapse.org/#!Synapse:syn53461705">Download the full aggregated counts matrix, rowData, and colData here</a>
</details>

1. **Get marker genes**:
    - Run `./04-get_marker_genes.ipynb`.

2. **Run celltype quality control and annotation**:
    - Run `./05-single_cell_qc_anno.ipynb` (Part I).
    - Run `./06-umaps.ipynb` (Part II).

3. **Save single cell data as SingleCellExperiment object**:
    - Run `./07-make_sce.ipynb`.

#### To Get K/L Gene Clusters:

<details>
<summary>Input Data</summary>
<a href="#">Download all the necessary input data through figshare (coming soon)</a>
</details>

1. **Run METIS and K/L algorithms**:
    - Run `./08-run_partitioning.py`.

2. **Benchmark clustering and partitioning methods**:
    - Run `./08-benchmarking_graph_partitioning.ipynb`.

#### To Perform Statistical Analyses:

<details>
<summary>Input Data</summary>
<a href="https://www.synapse.org/#!Synapse:syn53461705">Download the full annotated and QCed counts matrix, rowData, and colData here</a>

[Or download redacted versions (censored patient metadata) through the UCSC Single Cell Browser (coming soon)]
</details>

1. **Format data for input to stats analysis**:
    - Run `./09-stats_inputs.ipynb`.

2. **Compute gene scores and pathway enrichments**:
    - Run `./10-compute_stats.ipynb`.

3. **Gene score dimensionality reduction and clustering**:
    - Run `./11-projections.ipynb`.

4. **Format some data for plotting**:
    - Run `./13-plotting_inputs.ipynb`.

5. **Compute DEGs for the common ABCA7 variant**:
    - Run `./19-common_variant_analysis.ipynb`.

#### To Reproduce the Common Variant Analysis:

- Run `./20-common_var_plotting.ipynb` to plot the common variant analysis.

#### To Reproduce Lipidomic Results:

- Run `./16-lipidomics_PM.ipynb` to plot lipidomics aggregate data for the postmortem brain.

#### To Perform Additional Visualizations:

- Run `./12-KL_clusters.ipynb` to visualize graph partitioning results.
- Run `./14-figures.ipynb` to plot main figure panels.
- Run `./15-extended-figures.ipynb` to plot extended figures.
-

 Run `./17-variant_carrier_pie_charts.ipynb` to plot variant carrier proportions.
- Run `./18-specific_pathway_analysis.ipynb` to plot genes and pathways for targeted pathway analysis.
- Run `./21-basic_pie_charts.ipynb` to plot the common variant analysis.
- Run `./22-beta_ox_genes.ipynb` to plot the common variant analysis.

#### iPSC-Neuron Related

#### For All iPSC Neuronal Omics Analyses:

<details>
<summary>Input Data</summary>
<a href="#">Download all the necessary input data through figshare (coming soon)</a>
</details>

- Run `./23-seahorse.ipynb` to compute and plot oxygen consumption rates.
- Run `./24-lipidomics-iN.ipynb` to plot lipidomic related analysis results.
- Run `./25-metabolomics-iN.ipynb` to plot metabolome related analysis results.

#### For iPSC Neuronal Image Analyses:

- Visit the [confocalQuant GitHub Repository](https://github.com/djunamay/confocalQuant).

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
