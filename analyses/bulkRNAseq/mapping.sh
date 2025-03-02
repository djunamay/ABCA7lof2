#!/bin/bash
#SBATCH --job-name=make_index   
#SBATCH --time=4:00:00              
#SBATCH -c 40
#SBATCH --output=./bash_out/%x_%A_%a.out
#SBATCH --error=./bash_out/%x_%A_%a.err

#conda activate trim_env

for file in /home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/trimmed_fastq/*_1_sequence_val_1.fq; do
    base=$(basename ${file} 1_sequence_val_1.fq)
    STAR --genomeDir /home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/genome \
    --genomeLoad LoadAndKeep \
    --runThreadN 40 \
    --readFilesIn ${file} ${file/_1_sequence_val_1.fq/_2_sequence_val_2.fq} \
    --outFileNamePrefix /home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/mapped/${base}
done
