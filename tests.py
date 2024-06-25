import matplotlib.backends.backend_pdf
import numpy as np
from assertpy import assert_that
from ABCA7lof2.qc import compute_mito_fraction, get_fraction_mito, cluster_on_mito_and_counts, filter_cells, filter_on_gaussian_logliklihood
from ABCA7lof2.annotation import get_marker_indices, get_marker_matrix_norm, get_celltype_summary_score_per_cell, get_celltype_annotations_per_group
from numba_progress import ProgressBar
from datetime import datetime
import os

def test_get_fraction_mito():
    """
    Test function that returns fraction mito
    """
    gene_names = np.array(['Test1', 'MT-Test2', 'Test3', 'MT-Test4'])
    data = np.array(([4, 2, 3, 0], [0, 1, 2, 1]), dtype=int)
    mito_fractions, mito_index, total_counts = get_fraction_mito(gene_names, data, 'MT-')
    assert_that(np.array_equal(mito_fractions, np.array([2 / 9, 1 / 2])),
                f"Did not mito fractions correctly.").is_true()
    assert_that(np.array_equal(total_counts, np.array([9, 4])), f"Did not compute total counts correctly.").is_true()
    assert_that(np.array_equal(mito_index, np.array([False, True, False, True])),
                f"Did not get mito index correctly.").is_true()


def test_compute_mito_fraction():
    """
    Test computation of mito fraction and total counts
    """
    data = np.array(([0, 1, 0], [1, 0, 0], [0, 1, 2]), dtype=int)
    gene_index = [2]
    total_counts = np.empty(data.shape[0], dtype=int)
    mito_fractions = np.empty(data.shape[0], dtype=float)
    compute_mito_fraction(data, gene_index, mito_fractions, total_counts)
    ground_truth = np.array([0, 0, 2 / 3])
    assert_that(np.array_equal(mito_fractions, ground_truth), f"Did not compute mito fractions correctly.").is_true()
    ground_truth = np.array([1, 1, 3])
    assert_that(np.array_equal(total_counts, ground_truth), f"Did not compute total counts correctly.").is_true()


def test_cluster_on_mito_and_counts():
    """
    Test function that returns cluster assignments based on fraction mito and counts
    """
    N = 1000
    mito_fractions = np.random.uniform(size=N)
    total_counts = np.random.randint(0, 10000, size=N)
    result = cluster_on_mito_and_counts(mito_fractions, total_counts, sample_size=0.5)
    assert_that(result.shape).is_equal_to((N,))
    assert_that(result.mean()).is_greater_than_or_equal_to(0.5)


def test_filter_cells():
    """
    Test function that returns fraction mito
    """
    now=datetime.now().isoformat()
    path_to_outputs = './' + now + '_out'
    os.mkdir(path_to_outputs)

    N = 1000
    gene_names = np.array(['Test1', 'MT-Test2', 'Test3', 'MT-Test4'])
    counts2 = np.random.randint(low=0, high=1000, size=(N, len(gene_names)))
    meta = np.random.randint(low=0, high=1000, size=(N, 4))
    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/output.pdf")
    lower = 5
    upper = 2000
    mito_index, total_counts, mito_fractions, filtered_meta, filtered_counts = filter_cells(gene_names, counts2, meta,
                                                                                            0.1, lower, upper, pdf,
                                                                                            'MT-', path_to_outputs)
    assert_that(sum(total_counts >= lower)).is_equal_to(len(total_counts))
    assert_that(sum(total_counts <= upper)).is_equal_to(len(total_counts))

def test_get_marker_indices():
    marker_path = './test_data/adult_brain_markers.csv'
    gene_names = np.array(['MIR1302-2HG', 'ABLIM1', 'OR4F5'])
    marker_indices, celltype_dic = get_marker_indices(marker_path, gene_names)
    keys = set(celltype_dic.keys())
    assert_that(np.sum([x in keys for x in gene_names[marker_indices]])).is_equal_to(len(marker_indices))

def test_get_marker_matrix_norm():
    """
    Testing function that normalizes the matrix and returns only the marker gene matrix
    """
    filtered_counts = np.random.randint(low = 0, high = 50, size = (1000,10))
    marker_indices = np.random.randint(low = 0, high = 10, size = (1,3))[0]
    marker_out = np.memmap('marker_matrix_norm_test.npy', mode='w+', shape=(filtered_counts.shape[0], len(marker_indices)), dtype='float')
    total_counts = np.sum(filtered_counts, axis = 1)
    ground_truth = (filtered_counts/total_counts.reshape(-1,1))[:,marker_indices]
    with ProgressBar(total=filtered_counts.shape[0]) as numba_progress:
        get_marker_matrix_norm(filtered_counts, np.array(marker_indices), marker_out, total_counts, numba_progress)
    assert_that(np.array_equal(ground_truth, marker_out), f"Did not compute total marker matrix correctly.").is_true()

def test_get_celltype_summary_score_per_cell():
    """
    Testing function that scores celltype labels per-cell
    """
    unique_celltype_name = np.array(['cell_type_1', 'cell_type_2', 'cell_type_3'])
    celltype_name = np.array(['cell_type_1', 'cell_type_3', 'cell_type_3', 'cell_type_2', 'cell_type_1', 'cell_type_1', 'cell_type_3'])
    marker_out = np.random.rand(100, len(celltype_name))
    scores = get_celltype_summary_score_per_cell(unique_celltype_name, marker_out, celltype_name)
    assert_that(np.allclose(np.mean(marker_out[:,celltype_name=='cell_type_1'], axis = 1), scores[0]), f"Did not compute summary scores correctly.").is_true()

#def get_major_annotations():
    #unique_celltype_name = np.array(['cell_type_1', 'cell_type_2', 'cell_type_3'])
    #annotations_dic, summary_score = get_celltype_annotations_per_group(unique_celltype_name, celltype_summary_score, predictions)

def test_filter_on_gaussian_logliklihood():
    scores = np.random.rand(1000)
    out = filter_on_gaussian_logliklihood(scores)
    x = sum(out)
    y = len(out)-x
    assert_that(x>=y).is_true()


# write test for get_major_annotations
# check annotation correspondence to the manual annotation "ground truth"
# finish checking assign_major_celltypes functions
# where did filtering cells by gaussian confidence go? shouldn't that be before celltype specific QC --> or is it in the per-celltype QC?
# then do celltype-specific QC functions
# add intermediate QC plots
# plan the benchmarking --> what kinds of questions will give me confidence in the pipeline relative
# to other methods / speed generalizability across different tissues
# Questions: how many qc-ed cells (ground truth) make it in? how many unqced cells are left out (?)
# of the cells that commonly pass QC, what is the consensus in annotation between 'manual' and this pipeline?
# if the two points above show good accuracy --> is this method faster than other methods? what is the benefit to using it?

import numpy as np
from assertpy import assert_that
from ABCA7lof2.annotation import get_marker_matrix_norm, get_predictions, assign_major_celltypes


def test_get_marker_matrix_norm():
    """
    Testing function that normalizes the matrix
    """
    N = 1000
    marker_indices = np.array([0, 1, 0, 1])
    keep_cells = np.random.randint(0, 2, size=N)
    counts = np.random.randint(low=0, high=1000, size=(N, len(marker_indices)))
    total_counts = np.random.randint(10, 1000, size=N)
    marker_out = np.memmap('tests/marker_matrix_norm.npy', mode = 'w+', shape = (np.count_nonzero(keep_cells), np.count_nonzero(marker_indices)), dtype = 'float')
    get_marker_matrix_norm(counts, marker_indices, marker_out, total_counts, keep_cells)
    assert_that(marker_out.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(marker_out.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()

def test_get_predictions():
    N = 1000
    marker_indices = np.array([0, 1, 0, 1])
    keep_cells = np.random.randint(0, 2, size=N)
    counts = np.random.randint(low=0, high=1000, size=(N, len(marker_indices)))
    total_counts = np.random.randint(10, 1000, size=N)
    marker_out = np.memmap('tests/marker_matrix_norm.npy', mode = 'w+', shape = (np.count_nonzero(keep_cells), np.count_nonzero(marker_indices)), dtype = 'float')
    projected_matrix, predict, scores = get_predictions(counts, marker_indices, marker_out, total_counts, keep_cells, n_components_pca=2, sample_size=.5, n_components_gaussian=2)
    assert_that(projected_matrix.shape[1]).is_equal_to(np.count_nonzero(marker_indices)).is_true()
    assert_that(projected_matrix.shape[0]).is_equal_to(np.count_nonzero(keep_cells)).is_true()
    assert_that(len(predict)).is_equal_to(np.count_nonzero(keep_cells)).is_true()
    assert_that(len(scores)).is_equal_to(np.count_nonzero(keep_cells)).is_true()

def test_assign_major_celltypes():
    N = 1000
    marker_indices = np.array([1, 1, 1, 1, 0, 0])
    keep_cells = np.random.randint(0, 2, size=N)
    counts = np.random.randint(low=0, high=1000, size=(N, len(marker_indices)))
    total_counts = np.random.randint(10, 1000, size=N)
    marker_out = np.memmap('tests/marker_matrix_norm.npy', mode = 'w+', shape = (np.count_nonzero(keep_cells), np.count_nonzero(marker_indices)), dtype = 'float')
    projected_matrix, predict, scores = get_predictions(counts, marker_indices, marker_out, total_counts, keep_cells, n_components_pca=2, sample_size=.5, n_components_gaussian=2)
    marker_genes = ['M1', 'M2', 'M3', 'M4', 'M5']
    celltype_dic = dict(zip(marker_genes, ['C1', 'C1', 'C2', 'C2', 'C3']))
    marker_genes = ['M2', 'M3', 'M4', 'M1']
    annotations = assign_major_celltypes(marker_genes, celltype_dic, marker_out, predict)
    assert_that((len(scores))).is_equal_to(len(predict)).is_true()

def test_annotate_subtypes():
    N = len(annotation)
    annotation = np.array(['C1', 'C1', 'C2',
                                     'C1', 'C1', 'C2',
                                     'C1', 'C2', 'C2',
                                     'C1', 'C2', 'C2',
                                     'C1', 'C2', 'C2',
                                     'C1', 'C2', 'C2'])
    keep_cells = np.random.randint(0, 2, size=N)
    marker_indices = np.array([0, 1])
    counts = np.random.randint(low=0, high=1000, size=(N, 10))
    total_counts = np.random.randint(low=10, high=1000, size=N)
    n_components_pca = 2
    sample_size = 1
    n_components_gaussian = 2
    meta = np.random.randint(low=0, high=1000, size=(N, 5))
    meta_new = annotate_subtypes(keep_cells, annotation, meta, marker_indices, n_components_pca, n_components_gaussian, sample_size)
    # add assertions
    pass


import numpy as np
from assertpy import assert_that
from scqc2.qc import filter_genes, compute_mito_fraction, get_fraction_mito, cluster_on_mito_and_counts, filter_cells, filter_cells_by_celltype, filter_individuals, filter_genes_by_expression


def test_filter_genes():
    """
    Test that keeping genes with at least one nonzero entry
    """
    data = np.array(([0, 1, 0],[1, 0, 0]), dtype = int)
    keep_genes = np.full(data.shape[1], False, dtype=bool)
    filter_genes(data, keep_genes)
    ground_truth = np.array([True, True, False])
    assert_that(np.array_equal(keep_genes, ground_truth), f"Did not filter out genes with >0 counts in at least one sample.").is_true()


def test_compute_mito_fraction():
    """
    Test computation of mito fraction and total counts
    """
    data = np.array(([0, 1, 0],[1, 0, 0],[0, 1, 2]), dtype = int)
    gene_index = [2]
    total_counts = np.empty(data.shape[0], dtype=int)
    mito_fractions = np.empty(data.shape[0], dtype=float)
    compute_mito_fraction(data, gene_index, mito_fractions, total_counts)
    ground_truth = np.array([0, 0, 2/3])
    assert_that(np.array_equal(mito_fractions, ground_truth), f"Did not compute mito fractions correctly.").is_true()
    ground_truth = np.array([1, 1, 3])
    assert_that(np.array_equal(total_counts, ground_truth), f"Did not compute total counts correctly.").is_true()


def test_get_fraction_mito():
    """
    Test function that returns fraction mito
    """
    gene_names = np.array(['Test1', 'MT-Test2', 'Test3', 'MT-Test4'])
    data = np.array(([4, 2, 3, 0], [0, 1, 2, 1]), dtype=int)
    mito_fractions, mito_index, total_counts = get_fraction_mito(gene_names, data)
    assert_that(np.array_equal(mito_fractions, np.array([2/9, 1/2])), f"Did not mito fractions correctly.").is_true()
    assert_that(np.array_equal(total_counts, np.array([9, 4])), f"Did not compute total counts correctly.").is_true()
    assert_that(np.array_equal(mito_index, np.array([False, True, False, True])), f"Did not get mito index correctly.").is_true()


def test_cluster_on_mito_and_counts():
    """
    Test function that returns fraction mito
    """
    N=1000
    mito_fractions = np.random.uniform(size=N)
    total_counts = np.random.randint(0, 10000, size=N)
    result = cluster_on_mito_and_counts(mito_fractions, total_counts, sample_size=0.5)
    assert_that(result.shape).is_equal_to((N,))
    assert_that(result.mean()).is_greater_than_or_equal_to(0.5)


def filter_cells():
    """
    Test function that returns fraction mito
    """
    N=1000
    gene_names = np.array(['Test1', 'MT-Test2', 'Test3', 'MT-Test4'])
    counts = np.random.randint(low=0, high=1000, size=(N, len(gene_names)))
    keep_cells, mito_index, total_counts, mito_fractions = filter_cells(gene_names, counts, sample_size=0.1)
    assert_that(sum(total_counts[keep_cells])).is_greater_than_or_equal_to(1000*np.count_nonzero(keep_cells))
    assert_that(sum(total_counts[keep_cells])).is_less_than_or_equal_to(10000*np.count_nonzero(keep_cells))


def test_filter_cells_by_celltype():
    # filter cells by celltype
    mito_fractions = np.array([0.1, 0.1, 0.1, 0.1, 0.9, 0.9])
    total_counts = np.array([20, 20, 20, 20, 10, 10])
    sample_size = 1
    annotation = np.array(['C1', 'C2', 'C1', 'C2', 'C1', 'C2'])
    keep_cells_all = filter_cells_by_celltype(mito_fractions, total_counts, sample_size, annotation)
    ground_truth = np.array([1, 1, 1, 1, 0, 0]).astype('bool')

    # assertions
    assert_that(np.array_equal(keep_cells_all, ground_truth), f"Did not filter by celltype correctly.").is_true()
    mean_counts_keep = np.mean(total_counts[keep_cells_all])
    mean_counts_discard = np.mean(total_counts[np.logical_not(keep_cells_all)])
    assert_that(mean_counts_keep).is_greater_than(mean_counts_discard)
    mean_mito_keep = np.mean(mito_fractions[keep_cells_all])
    mean_mito_discard = np.mean(mito_fractions[np.logical_not(keep_cells_all)])
    assert_that(mean_mito_keep).is_less_than(mean_mito_discard)

def test_filter_individuals():
    celltype_annotations = np.array(['C1', 'C1', 'C2',
                                     'C1', 'C1', 'C2',
                                     'C1', 'C2', 'C2',
                                     'C1', 'C2', 'C2',
                                     'C1', 'C2', 'C2',
                                     'C1', 'C2', 'C2'])
    individual_annotation = np.array(['P1', 'P1', 'P1',
                                      'P2', 'P2', 'P2',
                                      'P3', 'P3', 'P3',
                                      'P4', 'P4', 'P4',
                                      'P5', 'P5', 'P5',
                                      'P6', 'P6', 'P6'])
    keep_cells = filter_individuals(celltype_annotations, individual_annotation)
    ground_truth = np.array(([0, 0, 0,
                            0, 0, 0,
                            1, 1, 1,
                            1, 1, 1,
                            1, 1, 1,
                            1, 1, 1])).astype('bool')
    assert_that(np.array_equal(keep_cells, ground_truth), f"Did not filter by celltype correctly.").is_true()


def test_filter_genes_by_expression():
    counts = np.array(([0, 1, 2, 3],
                      [0, 0, 0, 0],
                      [2, 2, 1, 0]))
    annotations = np.array(['C1', 'C1', 'C2', 'C2'])
    ground_truth = np.array(([1, 0, 1],
                            [1, 0, 1])).astype('bool')
    keep_genes = filter_genes_by_expression(annotations, counts)
    assert_that(np.array_equal(keep_genes, ground_truth), f"Did not filter genes by expression correctly.").is_true()
    pass
