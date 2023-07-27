#!/bin/bash

# SBATCH -J fastq_count
# SBATCH -t 0-48:00:00
# SBATCH --mail-type=ALL
# SBATCH --mail-user=djuna@mit.edu
# SBATCH --array=0-41
#SBATCH --mem=150G
#SBATCH -n 48
#SBATCH -N 1
# SBATCH -o /home/gridsan/djuna/github/ABCA7lof/raw_data/cellranger_counts_out/traces/slurm_ABCA7lof_counting_%A_%a.out
# SBATCH -e /home/gridsan/djuna/github/ABCA7lof/raw_data/cellranger_counts_out/traces/slurm_ABCA7lof_counting_%A_%a.err

# get mount path
mount1=$(/usr/local/bin/mount-squashfs.sh /home/gridsan/djuna/homer/github/archived_repos/ABCA7_LOF_2022/squash_files/batch_4819F.sqsh)
mount2=$(/usr/local/bin/mount-squashfs.sh /home/gridsan/djuna/homer/github/archived_repos/ABCA7_LOF_2022/squash_files/batch_4826F.sqsh)
mount3=$(/usr/local/bin/mount-squashfs.sh /home/gridsan/djuna/homer/github/archived_repos/ABCA7_LOF_2022/squash_files/171013Tsa.sqsh)

# run cellranger count
PATH=/home/gridsan/djuna/homer/github/archived_repos/scmod_R/f_cellranger/cellranger-6.1.2:$PATH 
eval $(sed -n "$SLURM_ARRAY_TASK_ID p" cellranger_count_executable.txt)
