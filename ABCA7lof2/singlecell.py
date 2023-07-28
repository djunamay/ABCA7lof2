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

#TODO: add docstrings to each function
def get_matrix_from_h5(filename):
    print('loading hd5..')
    #TODO: understand this function better - describe all the slots in the doc string
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
