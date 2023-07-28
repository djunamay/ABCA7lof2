import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import io
import scipy as scipy
import collections
import scipy.sparse as sp_sparse
import tables
import time
import numba as nb
from numba_progress import ProgressBar
import os
import ipdb 
import sys

def truncate_array(path, array, counter, Ncells, Ngenes):

    temp = open(path, 'r+')
    line = temp.readline()
    temp.close()

    constant = sys.getsizeof(line)

    temp = open(path, 'r+')
    
    if array.ndim==1:
        temp.truncate((array.itemsize*counter)+constant)
        line = temp.readline()
        temp.close()
        newline = line.replace(str(Ncells), str(counter))
    else:
        temp.truncate((array.itemsize*array.shape[1]*counter)+constant)
        line = temp.readline()
        temp.close()
        newline = line.replace(str((Ncells, Ngenes)), str((counter, Ngenes)))
        
    temp = open(path, 'r+')
    temp.writelines(newline)
    temp.close()

    
def get_matrix_from_h5(filename):
    '''
    function adapted from https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
    Args:
        filename (str) : path pointing to the counts matrix (genes x cell) output by the cellranger count function
        This file contains the counts for a single sample in a slot called 'matrix'
        Corresponding cell barcodes (column-names) are in the slot called 'barcodes' (must be unique)
    Returns:
        #TODO: specify what is being returned
    '''
    CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()
        return CountMatrix(feature_ref, barcodes, matrix)

def decode(numpy_array):
    return np.array([numpy_array[x].decode("utf-8") for x in range(len(numpy_array))])

def write_matrix(ID, counts_path_individual, cell_counts, cell_meta, subset_meta, counter, features_id_out, features_name_out, barcodes_out, Ngenes, i):
    #ipdb.set_trace()


    matrix = get_matrix_from_h5(counts_path_individual + ID + '/outs/filtered_feature_bc_matrix.h5')
    
    mat_array = matrix.matrix.T.toarray()
    temp_ncells = mat_array.shape[0]
    
    barcodes_out[counter:counter+temp_ncells]  = decode(np.array(matrix.barcodes))
    start = i*Ngenes
    stop = start+Ngenes
    features_id_out[start:stop] = decode(np.array(matrix.feature_ref['id']))
    features_name_out[start:stop] = decode(np.array(matrix.feature_ref['name']))
    
    write_counts_to_memmap(mat_array, cell_counts, counter)
    add_metadata(cell_meta, subset_meta, counter, temp_ncells)
    
    counter+=temp_ncells
    
    return counter

@nb.njit(parallel=True)
def write_counts_to_memmap(mat_array, cell_counts, counter):
    for x in nb.prange(mat_array.shape[0]):
        cell_counts[counter+x] = mat_array[x]

def add_metadata(cell_meta, subset_meta, counter, temp_ncells):
    for x in range(temp_ncells):
        cell_meta[counter+x] = subset_meta
        
def get_metadata(meta, Library_ID, barcodes):
    subset_meta = np.array(meta[meta['sample_id']==Library_ID])[0]
    subset_meta = np.array([subset_meta,]*len(barcodes))
    subset_meta = np.concatenate((barcodes[:,None], subset_meta), axis = 1)
    
    return subset_meta

def aggregate_fastqs(path_to_outputs, meta_path_individual, counts_path_individual, Ncells=300000, Ngenes=36601):
    print('preparing...')
    # initiate empty counts matrix and load as memmap
    counts_path = path_to_outputs + 'counts.npy'
    np.save(counts_path, np.empty(shape=(Ncells, Ngenes), dtype='int32'))
    cell_counts = np.lib.format.open_memmap(counts_path, mode='r+', dtype='int32', shape=(Ncells, Ngenes))

    # load metadata and get sample IDs
    ind_meta = pd.read_csv(meta_path_individual)
    n_metavars = ind_meta.shape[1]
    meta_path = path_to_outputs + 'metadata.npy'
    cell_meta = np.empty(shape=(Ncells, n_metavars), dtype=object)
    
    # initiate barcodes and features
    sample_ids = ind_meta['sample_id']
    n_ind = len(sample_ids)
    features_id_path = path_to_outputs + 'features_id.npy'
    np.save(features_id_path, np.empty(shape=Ngenes*n_ind, dtype='<U18'))
    features_id_out = np.lib.format.open_memmap(features_id_path, mode='r+', dtype='<U18')

    features_name_path = path_to_outputs + 'features_name.npy'
    np.save(features_name_path, np.empty(shape=Ngenes*n_ind, dtype='<U18'))
    features_name_out = np.lib.format.open_memmap(features_name_path, mode='r+', dtype='<U18')

    barcodes_path = path_to_outputs + 'barcodes.npy'
    np.save(barcodes_path, np.empty(shape=Ncells, dtype='<U18'))
    barcodes_out = np.lib.format.open_memmap(barcodes_path, mode='r+', dtype='<U18')
 
    # combine matrices
    counter = 0    
    counter = combine_matrices(sample_ids, cell_meta, ind_meta, cell_counts, counts_path_individual, features_id_out, features_name_out, barcodes_out, counter, Ngenes)
    
    # truncate arrays
    truncate_array(barcodes_path, barcodes_out, counter, Ncells, Ngenes)
    truncate_array(counts_path, cell_counts, counter, Ncells, Ngenes)

    # save metadata
    print('saving metadata')
    np.save(meta_path, cell_meta[:counter])
    print('done.')
    
def combine_matrices(sample_ids, cell_meta, ind_meta, cell_counts, counts_path_individual, features_id_out, features_name_out, barcodes_out, counter, Ngenes):
    
    for i in (range(len(sample_ids))):
        ID = sample_ids[i]
        print('aggregating ' + str(ID))
        subset_meta = np.array(ind_meta[ind_meta['sample_id']==ID])[0]
        counter = write_matrix(ID, counts_path_individual, cell_counts, cell_meta, subset_meta, counter, features_id_out, features_name_out, barcodes_out, Ngenes, i)
    return counter

