#!/bin/bash

# SBATCH -J sample_swap
# SBATCH --array=0-41
# SBATCH --mem=50
# SBATCH -o /home/gridsan/djuna/homer/github/ABCA7lof2/logs/slurm_ABCA7lof_counting_%A_%a.out
# SBATCH -e /home/gridsan/djuna/homer/github/ABCA7lof2/logs/slurm_ABCA7lof_counting_%A_%a.err

eval $(sed -n "$SLURM_ARRAY_TASK_ID p" sample_swap_exec.txt)
