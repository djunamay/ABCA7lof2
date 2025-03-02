#!/bin/bash
#SBATCH --job-name=make_index   
#SBATCH --time=4:00:00              
#SBATCH -c 40
#SBATCH --output=./bash_out/%x_%A_%a.out
#SBATCH --error=./bash_out/%x_%A_%a.err

#conda activate trim_env

featureCounts -a '/home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/reference/gencode.v47.primary_assembly.annotation.gtf' \
-o '/home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/featureCounts/counts.txt' \
-t 'exon' -g 'gene_id' \
--extraAttributes gene_name,gene_type \
-p \
-C \
-T 40 \
-B \ 
/home/gridsan/djuna/homer/github/ABCA7lof2/bulkRNAseq/mapped/*.sam

# -B: only count read pairs that have both ends aligned.
# -C: Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome but on different strands.
# -p: If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads; single-end reads are always counted as reads.
# -t: Specify the type of feature to count.