from argparse import ArgumentParser
import numpy as np
import numba as nb

#def get_test_data():

#@nb.njit(parallel=True) # why is numba not working here
def subset_mat(input_counts, output_counts, input_meta, output_meta, rand_index, progress_hook):
    for i in nb.prange(output_counts.shape[0]):
        output_counts[i] = input_counts[rand_index[i]]
        output_meta[i] = input_meta[rand_index[i]]
        progress_hook.update(1)

def get_data(meta_path, features_path, matrix_path): # must include the paths as arguments; also have to include them as user-specified paths to the aggregate argument
    meta = np.load(meta_path, allow_pickle=True)
    features = np.load(features_path, allow_pickle=True)
    counts = np.memmap(matrix_path, dtype='int16', shape=(meta.shape[0], features.shape[0]), mode='r')
    return meta, features, counts

def save_all(output_dir, filtered_meta_final, meta_names, keep_genes, gene_names, marker_indices, mito_index, projected_matrix, annotations): # what datatypes to save as?
    np.save(output_dir + '/filtered_metadata.npy', filtered_meta_final)
    np.save(output_dir + '/filtered_metadata_colnames.npy', meta_names)
    np.save(output_dir + '/keep_genes_per_celltype.npy', keep_genes)
    np.save(output_dir + '/gene_names.npy', gene_names)
    np.save(output_dir + '/marker_gene_index.npy', marker_indices)
    np.save(output_dir + '/mito_gene_index.npy', mito_index)
    np.save(output_dir + '/lower_dim_projection.npy', projected_matrix)
    if annotations.shape[0]>1:
       np.save(output_dir + '/inferred_marker_genes_per_cluster.npy', annotations)


def save_annotations(output_dir, filtered_meta_final, meta_names, projected_matrix, infer_markers):
    np.save(output_dir + '/metadata_with_annotations.npy', filtered_meta_final)
    np.save(output_dir + '/metadata_with_annotations_colnames.npy', meta_names)
    #np.save(output_dir + '/gene_names.npy', gene_names)
    #np.save(output_dir + '/marker_gene_index.npy', marker_indices)
    np.save(output_dir + '/lower_dim_projection.npy', projected_matrix)
    if infer_markers:
       np.save(output_dir + '/inferred_marker_genes_per_cluster.npy', annotations)


def save_data_by_celltype(meta, keep_genes, keep_cells, features):
    annotation = meta[:,-5]
    unique_celltype = np.unique(annotation)
    keep_cells = meta[:,-2]
    total_counts = meta[:,-1]
    for i in range(len(unique_celltype)):
        cellname = unique_celltype[i]
        keep_genes_temp = keep_genes[i]
        save_data(cellname, annotation, keep_cells, keep_genes_temp, total_counts, features, meta)

def save_data(annotation, meta, features, counts):
    index = annotation==cellname
    keep_cells_all = keep_cells & index
    norm_counts_out = np.memmap('normalized_counts_' + cellname + '.npy', mode='w+', shape=(len(keep_cells_all), len(keep_genes_temp)), dtype='float')
    counts_out = np.memmap('counts_' + cellname + '.npy', mode='w+', shape=(len(keep_cells_all), len(keep_genes_temp)), dtype='float')
    save_counts(norm_counts_out, counts_out, counts, keep_genes_temp, keep_cells_all, total_conts)
    np.save('features_' + cellname + '.npy', features[keep_genes_temp])
    np.save('features_' + cellname + '.npy', meta[keep_cells_all])

@nb.njit(parallel=True)
def save_counts(norm_counts_out, counts_out, counts, keep_genes, keep_cells_all, total_counts):
    for i in range(len(keep_cells_all)):
        counts_out[i] = counts[i][keep_genes] # make sure the mito genes are not included here
        norm_counts_out[i] = counts_out[i] / total_counts[i]

def process_metadata(scores, predict, projected_matrix, meta, mito_fractions, annotations, keep_cells, total_counts):
    meta = np.concatenate((meta, scores[:, None], predict[:, None], annotations[:, None], projected_matrix, mito_fractions[:, None], keep_cells[:, None], total_counts[:, None]), axis=1)
    return meta

def save_all_data(meta, keep_genes, counts, features):
    np.save('counts.npy', counts)
    np.save('features.npy', features)  # earlier, make sure that the matrices are the correct dimensions
    np.save('keep_genes.npy', keep_genes)
    np.save('meta.npy', meta)