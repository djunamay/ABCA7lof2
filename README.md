[![bioRxiv](https://img.shields.io/badge/bioRxiv-2024.07.28-b31b1b.svg?style=flat-square)](https://www.biorxiv.org/content/10.1101/2023.09.05.556135v2)

> **IMPORTANT**  
> This repo is under construction; Docstrings are still being added. Please check back soon!

# A Single-Cell Atlas of ABCA7 Loss-of-Function

This repository contains the main analysis code and links to raw and processed datasets, and relevant analysis pipelines, to reproduce (or extend on) results from our paper.

## Quickstart
- #### Clone the repository and install required packages:
    ```bash
    git clone git@github.com:djunamay/ABCA7lof2.git
    cd ABCA7lof2
    pip install -r requirements.txt
    ```
    - Ensure `mksquashfs` version 4.5 is installed on your system.

- #### See Data Availability below for links to the data.

- #### See Code Availability for a guide on reproducing the results. 
    - Each section is intended to be optionally run independently of the others. 
    - Information on accessing the input data for a given analysis is provided at the beginning of each section.

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

## ABCA7 loss-of-function snRNAseq QC & Annotation
> **<u>Data</u>**  
> - Follow [Get single cell BAM files](#get-single-cell-bam-files) above.
> - Then aggregate all the count files by running `./03-aggregate.ipynb`.
> - Alternatively, access the raw aggregated counts matrix, rowData, and colData [here](https://www.synapse.org/#!Synapse:syn53461705).
    
- Run `./04-get_marker_genes.ipynb` to get marker genes.
- Run `./05-single_cell_qc_anno.ipynb` to run celltype quality control and annotation.
- Run `./06-umaps.ipynb` to generate UMAPs.

## Compute Statistics
> **<u>Data</u>**  
> - Follow [ABCA7 loss-of-function snRNAseq QC & Annotation](#abca7-loss-of-function-snrna-qc--annotation) above.
> - Alternatively, access the QC'ed aggregated counts matrix, rowData, and colData [here](https://www.synapse.org/#!Synapse:syn53461705). 

- Run `./07-make_sce.ipynb` to save single cell data as SingleCellExperiment object.
- Run `./09-stats_inputs.ipynb` to format data for input to stats analysis.
- Run `./10-compute_stats.ipynb`.

## Gene and Pathway Statistics
> **<u>Data</u>**  
> - Follow [Compute Statistics](#compute-statistics) above.
> - Or, access the stats_input_data_0825.rds [here](https://figshare.com/s/c944697d9ec30ab06621).
- Run `./11-projections.ipynb`.
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

## Plot Results
> **<u>Data</u>**  
> Follow steps above or download figure input data.
- Run `./02-sample_swap.ipynb` to visualize sample swap results.
- Run `./14-figures.ipynb` to plot main figure panels.
- Run `./15-extended-figures.ipynb` to plot extended figures.    
- Run `./17-variant_carrier_pie_charts.ipynb` to plot variant carrier proportions.
- Run `./21-basic_pie_charts.ipynb` to plot the common variant analysis.

## Citation

Please cite our preprint using the following BibTeX entry if you use this code in your work:

```bibtex
@UNPUBLISHED{Von_Maydell2023-ce,
  title    = "Single-cell atlas of {ABCA7} loss-of-function reveals impaired
              neuronal respiration via choline-dependent lipid imbalances",
  author   = "von Maydell, Djuna and Wright, Shannon and Bonner, Julia Maeve
              and Staab, Colin and Spitaleri, Andrea and Liu, Liwang and Pao,
              Ping-Chieh and Yu, Chung Jong and Scannail, Aine Ni and Li,
              Mingpei and Boix, Carles A and Mathys, Hansruedi and Leclerc,
              Guillaume and Menchaca, Gloria Suella and Welch, Gwyneth and
              Graziosi, Agnese and Leary, Noelle and Samaan, George and Kellis,
              Manolis and Tsai, Li-Huei",
  abstract = "AbstractLoss-of-function (LoF) variants in the lipid transporter
              ABCA7 significantly increase the risk of Alzheimer's disease
              (odds ratio ∼2), yet the pathogenic mechanisms and the neural
              cell types affected by these variants remain largely unknown.
              Here, we performed single-nuclear RNA sequencing of 36
              humanpost-mortemsamples from the prefrontal cortex of 12 ABCA7
              LoF carriers and 24 matched non-carrier control individuals.
              ABCA7 LoF was associated with gene expression changes in all
              major cell types. Excitatory neurons, which expressed the highest
              levels of ABCA7, showed transcriptional changes related to lipid
              metabolism, mitochondrial function, cell cycle-related pathways,
              and synaptic signaling. ABCA7 LoF-associated transcriptional
              changes in neurons were similarly perturbed in carriers of the
              common AD missense variant ABCA7 p.Ala1527Gly (n = 240 controls,
              135 carriers), indicating that findings from our study may extend
              to large portions of the at-risk population. Consistent with
              ABCA7's function as a lipid exporter, lipidomic analysis of
              isogenic iPSC-derived neurons (iNs) revealed profound
              intracellular triglyceride accumulation in ABCA7 LoF, which was
              accompanied by a relative decrease in phosphatidylcholine
              abundance. Metabolomic and biochemical analyses of iNs further
              indicated that ABCA7 LoF was associated with disrupted
              mitochondrial bioenergetics that suggested impaired lipid
              breakdown by uncoupled respiration. Treatment of ABCA7 LoF iNs
              with CDP-choline (a rate-limiting precursor of
              phosphatidylcholine synthesis) reduced triglyceride accumulation
              and restored mitochondrial function, indicating that ABCA7
              LoF-induced phosphatidylcholine dyshomeostasis may directly
              disrupt mitochondrial metabolism of lipids. Treatment with
              CDP-choline also rescued intracellular amyloid$\beta$-42 levels
              in ABCA7 LoF iNs, further suggesting a link between ABCA7 LoF
              metabolic disruptions in neurons and AD pathology. This study
              provides a detailed transcriptomic atlas of ABCA7 LoF in the
              human brain and mechanistically links ABCA7 LoF-induced lipid
              perturbations to neuronal energy dyshomeostasis. In line with a
              growing body of evidence, our study highlights the central role
              of lipid metabolism in the etiology of Alzheimer's disease.",
  journal  = "bioRxiv",
  month    =  sep,
  year     =  2023
}

```
