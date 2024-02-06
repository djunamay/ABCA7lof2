#!/bin/bash

# get mount path
mount1=$(/usr/local/bin/mount-squashfs.sh /home/gridsan/djuna/homer/github/archived_repos/ABCA7_LOF_2022/squash_files/batch_4819F.sqsh)

# Define the source and destination directories
source_directory=$mount1/10x-fastqs-L2
destination_directory="/home/gridsan/djuna/homer/github/ABCA7lof2/open_data/synapse/fastq_files"

# Define the string to add to the filenames
string_to_add="_batch_4819F"

# Loop through all .gz files in the source directory
for file in "$source_directory"/*.gz; do
    # Extract the filename without the path
    filename=$(basename -- "$file")
    
    echo $filename
    
    # Since we're specifically dealing with .gz files, we know the extension already
    name="${filename%.fastq.gz}"

    # Construct the new filename with the added string
    new_filename="${name}${string_to_add}.fastq.gz"

    # Copy the file to the new location with the new filename
    cp "$file" "$destination_directory/$new_filename"
done
