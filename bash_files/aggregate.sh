#!/bin/bash

# SBATCH -J aggr
# SBATCH -t 0-48:00:00
# SBATCH -o /home/gridsan/djuna/homer/github/ABCA7lof2/logs/slurm_Jan2022_%A_%a.out
# SBATCH -e /home/gridsan/djuna/homer/github/ABCA7lof2/logs/slurm_Jan2022_%A_%a.err
# SBATCH --partition=xeon-g6-volta
# SBATCH --exclusive
# SBATCH --mem=280G

HDF5_USE_FILE_LOCKING='FALSE'

PATH=/home/gridsan/djuna/homer/github/archived_repos/scmod_R/f_cellranger/cellranger-6.1.2:$PATH 

cellranger aggr --id=ABCA7lof2 \
                --csv=/home/gridsan/djuna/homer/github/ABCA7lof2/raw_data/metadata/single_cell_individual_metadata.csv \
                --normalize=none \
                --nosecondary \
                --disable-ui \
                --jobmode=local \
                --localmem=280
