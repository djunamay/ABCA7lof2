## Software & Package Versions

### Core analysis tools
- **BMC/BCC pipelines** v1.8  
- **Cell Ranger (10x Genomics)** v6.1.2  
- **Harmony** (latest commit pulled 2025-06-23)

### Molecular-dynamics suite
- **PyMOL** v2.0  
- **CHARMM-GUI** (web service, release as of 2025-06-23)  
- **GROMACS** 2022.3  
- **VMD** v1.94  

### LC-MS data processing
- **LipidSearch** v4.2.27  
- **Compound Discoverer** v3.2  

### Metabolic-flux (Seahorse) analysis
- **XFe Assay software** v2.6.3.5  

### Electrophysiology & statistics
- **GraphPad Prism** v10  

### Runtime environments
- **Python** 3.8.13  
- **R** 4.2.3  

### Python packages
```text
scikit-learn==1.6.1      # GaussianMixture
statsmodels==0.14.4      # mixedlm
scanpy==1.10.3
umap-learn==0.5.7
gseapy==1.1.5
numpy==1.24.3            # Kernighan–Lin
networkx==3.2.1          # spring_layout
aicsimageio==4.14.0
````

### R packages

| Package | Version | Key functions used                    |
| ------- | ------- | ------------------------------------- |
| edgeR   | 4.4.0   | differential-expression               |
| limma   | 3.62.1  | `voom`, `lmFit`, `eBayes`, `topTable` |
| fgsea   | 1.32.2  | GSEA                                  |

### Additional bioinformatics utilities

* **METIS** ([https://github.com/KarypisLab/METIS](https://github.com/KarypisLab/METIS), commit 2025-06-23)
* **Trim Galore** v0.6.10
* **STAR aligner** v2.5.2b
* **featureCounts** v2.0.1

