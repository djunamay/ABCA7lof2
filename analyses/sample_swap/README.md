# Sample Swap Analysis 

This directory contains instructions and scripts for performing sample swap analysis on VCF files, specifically on VCF-BAM files from ROSMAP samples, and was used to perform sample-swap analysis in our paper.

## Data availability
- Download VCF files for chromosome(s) of interest following instructions here: [ROSMAPwgs GitHub](https://github.com/djunamay/ROSMAPwgs).
- We used `DEJ_11898_B01_GRM_WGS_2017-05-15_19.recalibrated_variants.vcf.gz.tbi` from `syn2580853` given ABCA7's location on Chr 19
- Generate BAM files for each sample: [snRNAseq_processing](../snRNAseq_processing/README.md)

### Genome Chain File (for coordinate conversion)

- **hg19ToHg38.over.chain.gz**  
  Source: [UCSC LiftOver hg19 → hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)

### Reference Genome (GRCh38)

- **GCF_000001405.26_GRCh38_genomic.fna.gz**  
  Source: [NCBI GRCh38 assembly (GCF_000001405.26)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26) 
## Requirements:
- **CrossMap (v0.6.5)**  
  [CrossMap Repository](https://github.com/Illumina/CrossMap)

- **HTSlib (v1.10.2)** for `bgzip` and `tabix`  
  [HTSlib Repository](https://github.com/samtools/htslib)

- **BCFtools (v1.10.2)**  
  [BCFtools Documentation](https://samtools.github.io/bcftools/bcftools.html)


## Code overview
- Follow these steps in order
- Adjust paths (`*/htslib-1.10.2/` and `*/bcftools`) according to your installation locations.
- Adjust data paths as needed in the respective files

### Remap the VCF file
```bash
./crossmap.sh
```
### Compress the remapped VCF
```bash
*/htslib-1.10.2/bgzip out.hg38.vcf --threads 20
```
### Sort compressed VCF
```bash
*/bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz
```
### Index sorted VCF
```bash
*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz
```
### Create `chr_name_conv.txt` 
```bash
cat <<EOT > chr_name_conv.txt
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
23 chrX
24 chrY
25 chrXY
26 chrM
EOT
```
### Rename chromosomes
```bash
*/bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
```
### Index annotated VCF
```bash
*/htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz
```

### Generate sample swap command for each sample
- see `./sample_swap_make_exec.ipynb` 

### Run sample swap
```bash
./sample_swap.sh
```

### Plot / inspect results
- See `./sample_swap.ipynb` notebook

## Resources

- [ROSMAPwgs GitHub](https://github.com/djunamay/ROSMAPwgs)
- [CrossMap](https://github.com/Illumina/CrossMap)
- [HTSlib](https://github.com/samtools/htslib)
- [BCFtools](https://samtools.github.io/bcftools/bcftools.html)


