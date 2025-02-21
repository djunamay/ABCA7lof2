#!/bin/bash
#SBATCH --job-name=trim   
#SBATCH --time=4:00:00              
#SBATCH -c 40
#SBATCH --output=./bash_out/%x_%A_%a.out
#SBATCH --error=./bash_out/%x_%A_%a.err

#conda activate trim_env

for file in /home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/raw_data/241121Tsa/*/*_1_sequence.fastq; do
    echo ${file}
    trim_galore --nextera --stringency 3 --cores 40 --paired --output_dir trimmed_fastq/ ${file} ${file/_1_sequence.fastq/_2_sequence.fastq}
done
