#!/bin/bash
#SBATCH --job-name=make_index   
#SBATCH --time=4:00:00              
#SBATCH -c 40
#SBATCH --output=./bash_out/%x_%A_%a.out
#SBATCH --error=./bash_out/%x_%A_%a.err

#conda activate trim_env

STAR --runMode genomeGenerate \
--genomeDir /home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/genome \
--genomeFastaFiles reference/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile reference/gencode.v47.primary_assembly.annotation.gtf \
--sjdbOverhang 75 \
--runThreadN 40
