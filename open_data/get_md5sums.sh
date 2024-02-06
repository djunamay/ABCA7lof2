# compute md5sums for fastq files 
source_directory="/home/gridsan/djuna/homer/github/ABCA7lof2/open_data/synapse/fastq_files"
output_file="/home/gridsan/djuna/homer/github/ABCA7lof2/open_data/synapse/fastq_files/md5sums.txt"

rm -f "$output_file"

cd "$source_directory"

for file in *.gz; do
    echo "$file"
    md5sum "$file" >> "$output_file"
done