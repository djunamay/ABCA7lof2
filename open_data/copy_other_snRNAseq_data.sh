# copy qc'ed counts and metadata

destination_directory="/home/gridsan/djuna/homer/github/ABCA7lof2/open_data/synapse/expression_matrices"

echo "gene_names"
cp "/home/gridsan/djuna/homer/github/ABCA7lof2/processed_data/single_cell/rowData.csv" "$destination_directory/qc_gene_names.csv"
echo "metadata"
cp "/home/gridsan/djuna/homer/github/ABCA7lof2/processed_data/single_cell/colData.csv" "$destination_directory/qc_sample_metadata.csv"
echo "counts"
cp "/home/gridsan/djuna/homer/github/ABCA7lof2/processed_data/single_cell/counts.mtx" "$destination_directory/qc_counts.mtx"


