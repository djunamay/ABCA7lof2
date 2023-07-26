/home/gridsan/djuna/Tools/QTLtools_1.1_Ubuntu16.04_x86_64 mbv --bam /home/gridsan/djuna/github/ABCA7lof/raw_data/cellranger_counts_out/D20-7432/outs/possorted_genome_bam.bam --out /home/gridsan/djuna/github/ABCA7lof/processed_data/VCF_BAM_comparison/test.txt --vcf /home/gridsan/djuna/github/ABCA7lof/raw_data/ROSMAP_WGS/DEJ_11898_B01_GRM_WGS_2017-05-15_19.recalibrated_variants.vcf.gz --filter-mapping-quality 150 --reg chr19:1039997-1065572

/home/gridsan/djuna/Tools/QTLtools_1.1_Ubuntu16.04_x86_64 mbv --bam /home/gridsan/djuna/github/ABCA7lof/raw_data/cellranger_counts_out/D20-7432/outs/possorted_genome_bam.bam --out /home/gridsan/djuna/github/ABCA7lof/processed_data/VCF_BAM_comparison/test.txt --vcf /home/gridsan/djuna/github/ABCA7lof/raw_data/ROSMAP_WGS/DEJ_11898_B01_GRM_WGS_2017-05-15_19.recalibrated_variants.vcf.gz --reg 19:500000-1500000

# first use crossmap to remap the vcf file (see crossmap.sh)
# then compress with bgzip: /home/gridsan/djuna/Tools/htslib-1.10.2/bgzip out.hg38.vcf --threads 20
# sort the vcf file: /home/gridsan/djuna/genome_apps/bin/bcftools sort out.hg38.vcf.gz -o out.hg38.sorted.vcf.gz
# then generate the corresponding tabix file /home/gridsan/djuna/Tools/htslib-1.10.2/tabix -p vcf out.hg38.sorted.vcf.gz
# check how the chromosomes are named in the BAM file: /home/gridsan/djuna/genome_apps/bin/samtools idxstats /home/gridsan/djuna/github/ABCA7lof/raw_data/cellranger_counts_out/D20-7432/outs/possorted_genome_bam.bam | head -n 3
# /home/gridsan/djuna/genome_apps/bin/bcftools annotate --rename-chrs chr_name_conv.txt out.hg38.sorted.vcf.gz -Oz -o out.hg38.sorted.ChrNamed.vcf.gz --threads 40
# check if the chromosome names changed: /home/gridsan/djuna/genome_apps/bin/bcftools query -f '%CHROM %POS %REF %ALT\n' out.hg38.sorted.ChrNamed.vcf.gz | head -3
# then generate the corresponding tabix file /home/gridsan/djuna/Tools/htslib-1.10.2/tabix -p vcf out.hg38.sorted.ChrNamed.vcf.gz

### double-check that the original vcf file was definitly in Hg37 format
/home/gridsan/djuna/Tools/QTLtools_1.1_Ubuntu16.04_x86_64 mbv --bam /home/gridsan/djuna/github/ABCA7lof/raw_data/cellranger_counts_out/D20-7432/outs/possorted_genome_bam.bam --out /home/gridsan/djuna/github/ABCA7lof/processed_data/VCF_BAM_comparison/mbv_out_D207432.txt --vcf /home/gridsan/djuna/github/ABCA7lof/processed_data/VCF_BAM_comparison/out.hg38.sorted.ChrNamed.vcf.gz #--reg chr19:500000-1500000

