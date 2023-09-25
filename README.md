[![bioRxiv](https://img.shields.io/badge/bioRxiv-202023.09.05-b31b1b.svg?style=flat-square)](https://www.biorxiv.org/content/10.1101/2023.09.05.556135v1)

[This repo is under construction]
# A single-cell atlas of ABCA7 loss-of-function reveals lipid disruptions, mitochondrial dysfunction and DNA damage in neurons

In our [paper](https://www.biorxiv.org/content/10.1101/2023.09.05.556135v1), we explored the impact of ABCA7 loss-of-function (LoF) variants on Alzheimer's disease (AD) by conducting single-nuclear RNA sequencing on 36 human post-mortem samples from the prefrontal cortex. We found that ABCA7 LoF variants resulted in gene expression changes across all major cell types, with the most marked alterations observed in excitatory neurons. These changes influenced lipid metabolism, mitochondrial function, DNA damage, and NF-kB signaling. Through functional assays, we confirmed elevated levels of mitochondrial dysfunction, DNA damage, and NF-kB activation in neurons. We also used mass spectrometry to reveal that ABCA7 LoF led to significant changes in the lipidome, including increased triglycerides and altered phospholipid species. Our study provides a detailed transcriptional atlas of ABCA7 LoF effects in the human brain and suggests that lipid dysregulation in neurons may be a key mechanism by which ABCA7 LoF increases the risk of AD. This GitHub repository contains all the code necessary to replicate our analyses and to explore, generate, and analyze this detailed transcriptional atlas.

## Data

Follow these instructions to access the data generated and used as part of this study.

- For the processed **WGS data**, follow instructions in the [ROSMAPwgs](https://github.com/djunamay/ROSMAPwgs) repository
> Then run the main.py script with the following parameters
> ```bash
> python main.py --outdir './raw_data/ROSMAP_WGS' --username <USERNAME> --pw <PASSWORD> --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extension 'recalibrated_variants.vcf.gz' --extract_HIGHandMED_annotations False --download True
> python main.py --outdir './raw_data/ROSMAP_WGS' --username <USERNAME> --pw <PASSWORD> --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extension 'annotated.coding.txt' --extract_HIGHandMED_annotations False --download True
> python main.py --outdir './raw_data/ROSMAP_WGS' --gene_list "['SORL1', 'TREM2', 'ABCA7', 'ATP8B4', 'ABCA1', 'ADAM10']" --extract_HIGHandMED_annotations True --download False
> ```
- For the processed and raw **post-mortem snRNAseq data**, go to [Synapse](https://linktosynapse) to request the data

- For the raw **post-mortem lipidomic data**, go to [Synapse](https://linktosynapse) to request the data

- For **post-mortem lipidomic stats results** and **post-mortem snRNAseq stats results**, go to [OSF](https://osf.io/vn7w2/) to download

- For raw and processed **post-mortem imaging data**, go to [OSF](https://osf.io/vn7w2/) to download

- For all processed and raw **iPSC data**, go to [OSF](https://osf.io/vn7w2/) to download
  
## Analyses

### run cellranger counting
> <details>
> <summary>Show the methods</summary>
> Library demultiplexing was performed using the BMC/BCC pipelines (https://openwetware.org/wiki/BioMicroCenter:Software). Fast-q reads were aligned to human genome GRCh38 and counted using the cellranger count() function from cellranger version 6.1.2. (10x Genomics). Introns were included in counting, to allow for detection of unspliced transcripts and the expected number of cells was set to 5000. Otherwise cellranger (v.6.1.2) default parameters were used. Counts across individual samples were then aggregated using a custom aggregation script, resulting in a total of 150,456 cells. 
> </details>
>
> <details>
> <summary>1. Download the data</summary>
>
> [Download FASTQ files here](https://linktosynapse)    
> </details>
>
> <details>
> <summary>2. make the squash file system</summary>
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
> <summary>3. run cellranger counting</summary>
>
> ```bash
> # count the FASTQ files:
> sbatch --array 1-42 */bash_files/cellranger_count.sh
> */bash_files/check_success.sh # iterate over all logs and check whether pipeline was successful before moving to aggregation
> ```
> </details>
>
> 4. run *`./03-aggregate.ipynb`* to aggregate all the count files

### sample swap
> <details>
> <summary>Show the methods</summary>
> Sample swap analysis was performed using a previously established pipeline (MVV; QTLtools_1.1) (2), which compares allelic concordance between genomic and transcriptomic sequencing data. As input, we used the BAM files generated in the cellranger counting step and the chromosome 19 (the chromosome harboring ABCA7) variant call files (VCF). When comparing the concordance of BAM and VCF data for homozygous and heterozygous sites, the expected WGS sample should appear as a clear outlier.
> </details>
>
> <details>
> <summary>1. Download the data</summary>
>
> See sections **`run cellranger counting`**  and **`processed WGS data`** above to get BAM files and WGS data.
> </details>
>
> 2. run *`/bash_files/crossmap.sh`* to remap the vcf file
> <details>
> <summary>3.then run the following commands in bash:</summary>
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
> 4. run *`sample_swap_make_exec.ipynb`* to make text file to iterate through 
> 5. run *`sample_swap.sh`* to run sample swap 
> 6. run *`./02-sample_swap.ipynb`* to visualize sample swap results 

### celltype annotation & QC
> <details>
> <summary>Show the methods</summary>
> 
> **Cell filtering metrics**    
> * Prior to cell type annotation, we performed a series of quality control steps on the aggregated counts matrix. First, we removed cells, for which the number of detected genes (Ng) did not fall in the interval [500, 10000], where Ng is defined for each cell as the number of genes, where counts >0. 
> * Next, we removed all cells with a high fraction of counts coming from mitochondrial-encoded genes. Mitochondrial fraction (Mf) is a commonly used per-cell metric to measure compromised nuclear integrity, where high fractions indicate low-quality nuclei, where Mf = (total counts mitochondrially-encoded genes)/(total counts all genes). We fit a gaussian mixture model to log2-transformed per-cell Mf values, using grid search to choose the optimal number of components and covariance type. The model with five components and full covariance had the lowest bayesian information criterion (BIC) score. Cells belonging to the component with the highest average Mf scores were presumed to constitute a population of low-quality cells and removed from further analysis. This initial filtering removed circa 20,000 cells.   
> * Considering all remaining cells in marker-gene expression space, where marker genes include only known cell type-specific genes for the major human PFC cell type, including astrocytes (with 159 markers), excitatory neurons (113 markers), inhibitory neurons (83 markers), microglia (97 markers),  oligodendrocytes (179 markers), OPCs (143 markers), and  vascular cells (124 markers) (Reference 1; Table S2) normalized to total library counts, we performed incremental PCA (IncrementalPCA from sklearn.decomposition) on this mean-centered standardized matrix to project cells from the marker gene space onto the top 50 principal components sorted by variance. Visually, the cells formed a number of gaussian-like clusters when the first two principal components were examined. Under the assumption that each gaussian cluster represented a different cell type in the brain, we next fit a gaussian mixture model (GaussianMixture from sklearn.mixture) to the projected data, using grid search (GridSearchCV from sklearn.preprocessing) to choose the optimal number of components and covariance type. The model with ten components and full covariance had the lowest BIC score. Indeed, each resulting cell cluster was robustly enriched for a subset of major cell type markers in the brain, indicating a cluster of astrocytes, microglia, OPCs, oligodendrocytes, excitatory neurons, and inhibitory neurons, and a heterogeneous cluster of vascular cells. 
> * To remove cells that were not well-captured by this model and likely represent low-quality cells, we next computed the per-cell logliklihood (i.e. the liklihood of the observed data, given the model) and removed cells with a liklihood \< -100. We also removed two gaussian clusters whose liklihood distributions constituted clear outliers compared to remaining clusters. The excluded cells had significantly lower total counts and higher mitochondrial fractions compared to those that passed the liklihood filter, suggesting that the removed cells indeed were low quality. When examining the data visually projected onto the first two principal components, this filtering removed many of the cells that were not visibly associated with a main gaussian cluster. Together, this filtering removed an additional circa 12,000 cells, leaving a total of 118,668 cells. 
>
> **Gene filtering metrics**
> * For the remaining downstream analysis we only considered genes that were both nuclear-encoded and protein-coding, which constituted a total of 19384 genes, based on annotation of ensembl GRCh38p12. 
>
> **Cell type annotations**
> * To avoid biased cell type annotations due to technical artifacts associated with sequencing batch and individual-of origin, we first applied the Python implementation of the Harmony algorithm (3) with individual-of-origin as indicator vector to the low-dimensional embedding of cells (first 50 principal components) remaining after the initial rounds of quality control (see above Cell filtering metrics). Next, we computed a neighborhood graph on the Harmony-corrected values in the PC embedding space, as implemented in the Scanpy (4) Python package, using default parameters. Finally, we applied the Leiden graph-clustering algorithm to cluster this neighborhood graph of cells, using the Scanpy implementation of the Leiden algorithm (5).
> * We used the Scanpy ‘rank_genes_groups’ function to compute top marker genes per Leiden cluster. Internally, this function uses a T-test to compute the relative enrichments of genes for each Leiden cluster compared to all other Leiden clusters. We assigned a major cell type label (‘Ex’, ‘In’, ‘Ast’, ‘Mic’, ‘Oli’, or ‘Opc’, ‘Vascular’) by computing per-cluster average log2-fold-changes (logFC) for respective cell type markers (Reference 1; Table S2) and assigning the cell type with the highest logFC, where large and positive logFCs indicate high relative expression of a gene in a given Leiden clusters compared to all other clusters. 
> * Finally, we sub-clustered cells from each major cell type using the Leiden clustering algorithm and examined distributions of mitochondrial fractions and total counts among subclusters of the same cell type. Clusters whose mean mitochondrial fraction was >2 standard deviations (sd) above the mean or whose mean total counts were < 2 sd below the mean or >2 sd above the mean (when comparing sub-clusters of a single cell type) were removed. Manual inspection of the removed clusters revealed that they tended to have fewer cells and low individual-level representations, and were not well-connected in the graph.
>
> **Individual-level filtering** 
> * After all rounds of qc as described above, we noted a subset of individuals (N=6) with very few cells (<500) and these subjects were removed from further analysis, resulting in 24 control individuals and 12 ABCA7 LoF individuals. None of these individuals carried ABCA7 PTC variants and removing them did not substantially alter metadata distributions. 
> </details>
> <details>
> <summary>1. Download the data</summary>
>
> [Download the aggregated counts matrix, rowData, and colData here](https://linktosynapse)    
> </details>
>
> 2. run *`./04-get_marker_genes.ipynb`* to get marker genes for celltype annotation 
> 3. run *`./05-single_cell_qc_anno.ipynb`* to run celltype quality control and annotation I 
> 4. run *`./06-umaps.ipynb`* to run celltype quality control and annotation II 
> 5. run *`./07-make_sce.ipynb`* to save single cell data as singlecellexperiment object 
    
### gene clusters 
> <details>
> <summary>Show the methods</summary>
>
> * To reduce the solution's computational search space, we reformulated the gene-pathway association problem as a bipartite graph G constructed from all the genes in the Leading Edge subset (LE) and their associated pathways. LE was defined as the set of 268 genes driving the enrichment signal for pathways that passed a significance threshold of p<0.05 (FGSEA) in Con vs ABCA7 LoF excitatory neurons. G was constructed from an n x m unweighted adjacency matrix, where n represented the number of LE genes and m the number of pathways associated with four or more LE genes, as specified in the WikiPathway database.
> * We chose to group gene-pathways into clusters of approximately equal size, making this a graph partitioning problem, because we found that removing this constraint made the grouping results highly susceptible to outliers (Supplementary Text; Fig. S8C). Of the three graph partitioning algorithms tried, METIS and the Kernighan-Lin (K/L) algorithms had the lowest loss (Supplementary Text; Fig. S8B). Both METIS and K/L achieved very comparable losses (within 1.8% of each other, after 5.0x10e4 random initiations) and came to almost identical solutions (Rand index=0.98, after 5.0x10e4 random initiations) (Supplementary Text; Fig. S8B,D-F). We proceeded with the K/L algorithm for gene-pathway groupings as we found this algorithm to perform consistently better than METIS across a wider range of graph sizes (not shown).
> * The K/L algorithm was implemented in Python based on its original paper (14) and run with parameters set as C=0, KL_modified=True, random_labels=True, unweighted=True, and K=50 to partition G into 8 groups. We performed 5.0x10e4 random initiations on G and report the partitioning with the lowest loss among all initiations.
> * For benchmarking results, see the correpsonding [github repo](https://github.com/djunamay/geneclusters).
> </details>
> <details>
> <summary>1. Get the data</summary>
>
> [Download the gene-pathway matrix here](https://osf.io/vn7w2/)    
> </details>
>
> 1. run *`./08-run_partitioning.py`* to run METIS and K\L algorithms 
> 2. see *`./08-processing_gsets.ipynb`* to benchmark clustering and partitioning methods

### stats
> <details>
> <summary>Show the methods</summary>
>
> **Differential gene expression**
> * Per-gene count values were summed for each cell-type-by-individual combination, resulting in 36 pseudobulk gene expression vectors for each of the six major cell types. For each cell type, only genes with a nonzero detection rate >0.10 were considered for differential expression. Pseudobulk counts were normalized using the edgeR (6, 7) TMM method. The residual mean-variance trend not explained by the multivariate linear model (formalized below), was removed using Limma-Voom (8). Unknown sources of variance were captured in the model using surrogate variable analysis (SVA) (9). Limma’s lmfit, eBayes, and toptable functions (10) were then used to estimate differential gene expression statistics, as reported in Data S3. 
> * The following model was fit for each cell type:
> * Gi = 𝛃0 * ABCA7LoF + 𝛃1 * msex + 𝛃2 * nft + 𝛃3 * amyloid + 𝛃4 * age_death + 𝛃5*PMI + 𝛃6 * batch + 𝛃7 * APOE4 +  𝛃8 * SV0
> * Gi refers to a vector of expression profiles of size 1 x 36 for a gene i in a given cell type. ABCA7LoF is a binary variable, encoding the presence of an ABCA7 variant predicted to cause LoF (see Data S1). See Supplementary Text for descriptions of the remaining variables included in the model. SV0 refers to the first surrogate variable estimated from the data. The exact number of surrogate variables per cell type to include as additive terms in the model was estimated using the num.sv() function in R. Up to 10 SVs were included.
> 
> **Gene-set enrichment analysis**
> * Genes were rank ordered based on their scores S (see description in Gene-pathway projections). An R implementation of gene set enrichment analysis (9, 13) (fast gene set enrichment analysis, fGSEA (11)) was run with 10,000 permutations to estimate the statistical overrepresentation of gene sets in the WikiPathways database (Table S2) within high-scoring (|S|), differentially expressed genes. Gene sets with a minimum size of 5 and a maximum size of 1000 were considered. 
> * To query a more comprehensive set of specific biological themes related to lipids and NF-kB, we ran fGSEA on a larger database of pathways (including HumanCyc, KEGG, Reactome, Biocarta, WikiPathways, and GO BP), after filtering these pathway databases for lipid ('sterol', 'lipid', 'glycer', 'fatt', 'ceramide', 'phosphatidyl') and NFkB (‘kappa’) -related terms, respectively.	
> 
> **Overlap with CRISPRi perturb-sequencing dataset**
> * Glutamatergic Neuron-RNA-Seq-CRISPRi (2020) data differential gene expression statistics (from CRISPRi perturbed vs non-perturbed cells) were downloaded from the online crisprbrain.org resource (Table 2). For each CRISPRi gene target, downstream gene expression changes were summarized as scores (where Score=sign(log2(FC))*-log10(p-value)), for all genes, for which average expression log2CPM > 0. FGSEA was used (parameters: minSize = 5, maxSize = 1000, nPermSimple=10000) to compute the enrichment of all excitatory neuronal K\L gene clusters (Data S7) among highly differentially expressed genes (Score) per CRISPRi target. Only CRISPRi target genes that had a nonzero detection rate >0.25 in post-mortem excitatory neurons, as assessed by snRNAseq, and had an ABCA7 LoF perturbation score |S| > 0.5 in excitatory neurons, as assessed by snRNAseq in the post-mortem human brain, were considered for this analysis. 
>
> **Lipidomic data analysis**
> * Lipids were identified and their signal integrated using the Lipidsearch © software (version 4.2.27, Mitsui Knowledge Industry, University of Tokyo). Integrations and peak quality were curated manually before exporting. Statistical significance for differences in peak distributions between control and ABCA7 LoF or WT and p.Glu50fs*3 were computed by two-sided unpaired T-test. 
> 
> **Gene-pathway projections**
> * For each cell type, we computed a set of gene-wise scores quantifying the direction and statistical significance of gene expression changes (computed as part of the differential gene expression analysis) associated with ABCA7 LoF: S = sign(log2FC) * -log10(p-value), where log2FC>0 indicates up-regulation in ABCA7 LoF vs control. Top differentially expressed genes per cell type (|S|>1.3) were projected from 6-dimensional score space, where each dimension captures ABCA7 LoF perturbation scores in one of the major cell types (Ex, In, Ast, Mic, Oli, OPC), into two dimensions, using the UMAP algorithm (11, 12) (using the umap Python package). Gene scores that were not detected >10% in a given cell type were set to 0. We performed a grid search for gaussian mixture parameters (parameter 1: N components; parameter 2: covariance type) on the embedded cells (using the Python sklearn package) to assign genes to clusters in the 2D space. We proceeded with the model with the lowest BIC score, which had 15 components and a tied covariance matrix. 
> * Each cluster was assigned representative pathway names by testing genes in that cluster for enrichment with Gene Ontology Biological Process pathways (Table S2) against the background of all genes in the embedding space, by hyper-geometric enrichment (using the Python package gseapy). Pathways with an enrichment p-value < 0.01 were considered for cluster annotation. 
Per-cell-type perturbation scores (Sc) for each cluster were computed as the average gene score S (for a given cell type) for all genes in that cluster. The statistical significance of each cell type-specific cluster score was assessed by permuting cluster assignments (100,000 permutations). 
>
> </details>
> <details>
> <summary>1. Get the data</summary>
> </details>
>
> 2. see *`./09-stats_inputs.ipynb`* to format data for input to stats analysis 
> 3. see *`./10-compute_stats.ipynb`* to compute gene scores and pathway enrichments 
> 4. see *`./11-projections.ipynb`* for gene score dimensionality reduction and clustering 
> 5. see *`./13-plotting_inputs.ipynb`* to format some data for plotting 

### plots
> <details>
> <summary>Show the methods</summary>
> 
> **gene-pathway graphs**
> * Gene-pathway graph layouts were computed using the networkx Python package using the spring layout algorithm, with 10,000 iterations. Layouts were visualized using the matplotlib pyplot package in Python. 
> * Representative pathways for each cluster were inferred from the graph, by averaging the ABCA7 LoF perturbation scores S for all genes in the cluster of interest sharing an edge with the pathway in question. Scores for pathways with intra-cluster degrees>=5 were reported in the figures. Manually picked subsets of genes with the largest scores (|S|>1) were reported in the figures. All gene statistics are reported in Data S3 and cluster assignments are reported in Data S7. 
> </details>
> <details>
> <summary>1. Get the data</summary>
> </details>  
>
> 2. see *`./12-KL_clusters.ipynb`* to visualize graph partitioning results 
> 3. see *`./14-figures.ipynb`* to plot main figure panels 
> 4. see *`./15-extended-figures.ipynb`* to plot extended figures

## Citation
If you use this code in your work, please cite using the following BibTeX entry:
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
