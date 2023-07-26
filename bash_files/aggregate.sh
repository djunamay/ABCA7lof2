#!/bin/bash

# SBATCH -J fastq_count
# SBATCH -t 0-48:00:00
# SBATCH -o /home/gridsan/djuna/github/ABCA7_LOF_2022/pre-processing/cellranger_count/cellranger_count_output/out/slurm_Jan2022_%A_%a.out
# SBATCH -e /home/gridsan/djuna/github/ABCA7_LOF_2022/pre-processing/cellranger_count/cellranger_count_output/err/slurm_Jan2022_%A_%a.err
# SBATCH --mail-type=ALL
# SBATCH --mail-user=djuna@mit.edu
# SBATCH --mem=100G
# SBATCH --cpus-per-task=40
# SBATCH --exclusive 

PATH=/home/gridsan/djuna/github/scmod_R/f_cellranger/cellranger-6.1.2:$PATH 

cellranger aggr --id=ABCA7_LOF_Jan12_2022 \
                --csv=metadata_for_aggregation.csv \
                --normalize=none \
                --nosecondary \
                --localcores=40 \
                --disable-ui \
                --jobmode=local \
                --localmem=50
