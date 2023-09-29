import numba as nb
import numpy as np
import pandas as pd
from sklearn.decomposition import IncrementalPCA
from sklearn.mixture import GaussianMixture
from ABCA7lof2.qc import filter_cells_by_major_annotation
from ABCA7lof2.setup import process_metadata
from numba_progress import ProgressBar
import localreg
from ABCA7lof2.qc import gmm_bic_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler

def get_marker_indices(marker_path, gene_names):
    """
    Retrieve indices of celltype-marker genes and create a celltype-marker gene dictionary.

    Args:
        marker_path (str):
            Path to a CSV file containing marker gene annotations with the following columns:
            - marker: Gene names (must intersect with the input gene_names argument).
            - major_celltype: Corresponding cell type.

        gene_names (numpy.ndarray):
            1D array containing gene names with '<U17' data type; array length = N-features.

    Returns:
        tuple:
            A tuple containing two elements:
            - marker_indices (list): List of indices for celltype-marker genes in the gene_names array.
            - celltype_dic (dict): Dictionary mapping marker genes to their corresponding cell types.

    Notes:
    ------
    This function reads marker gene annotations from a CSV file, retrieves the indices of celltype-marker
    genes in the input gene_names array, and creates a dictionary mapping marker genes to cell types.

    Example:
    --------
    >>> marker_path = "marker_annotations.csv"
    >>> gene_names = np.array(['GeneA', 'GeneB', 'GeneC'], dtype='<U17')
    >>> marker_indices, celltype_dic = get_marker_indices(marker_path, gene_names)
    >>> marker_indices
    [0, 1]  # Indices of celltype-marker genes in the gene_names array.
    >>> celltype_dic
    {'GeneA': 'CellTypeA', 'GeneB': 'CellTypeB'}  # Dictionary mapping marker genes to cell types.
    """
    print('getting marker indices...')
    markers = pd.read_csv(marker_path)
    dic = dict(zip(gene_names, range(0, len(gene_names))))
    targets = np.unique(np.array(markers['marker']))
    targets = targets[[targets[x] in set(gene_names) for x in range(len(targets))]]
    marker_indices = [dic[targets[x]] for x in range(len(targets))]
    celltype_dic = dict(zip(markers['marker'], markers['major_celltype']))
    return marker_indices, celltype_dic

@nb.njit(parallel=True)
def get_marker_matrix_norm(filtered_counts, marker_indices, marker_out, total_counts, progress_hook):
    """
    Calculate and return the normalized matrix of marker genes for cell-type annotation.

    Args:
        filtered_counts (numpy.ndarray):
            2D array of counts post-cell filtering (N-cells x N-features).
        marker_indices (list):
            List of column indices indicating marker genes in the filtered_counts matrix.
        marker_out (numpy.memmap):
            Memory-mapped array to store the normalized marker matrix (N-cells x len(marker_indices)).
        total_counts (numpy.ndarray):
            1D array of per-cell total counts, length = N-cells.
        progress_hook (ProgressHook):
            An object for tracking progress during computation.

    Returns:
        None

    Notes:
    ------
    This function calculates the normalized matrix of marker genes for cell-type annotation.
    It divides the counts of marker genes by the total counts for each cell and stores
    the result in a memory-mapped array.

    Example:
    --------
    >>> filtered_counts = np.array([[10, 5, 3], [8, 4, 2]], dtype='int32')
    >>> marker_indices = [0, 2]
    >>> marker_out = np.memmap("marker_normalized.npy", dtype='float32', mode='w+', shape=(2, 2))
    >>> total_counts = np.array([18, 14], dtype='int32')
    >>> progress_hook = ProgressHook(total=2)
    >>> get_marker_matrix_norm(filtered_counts, marker_indices, marker_out, total_counts, progress_hook)
    # Calculates and stores the normalized marker matrix.
    """
    for i in nb.prange(filtered_counts.shape[0]):
        marker_out[i] = (filtered_counts[i][(marker_indices)])/total_counts[i]
        progress_hook.update(1)

def run_ipca(marker_mat, n_components_pca, sample_size):
    """
    Run Incremental Principal Component Analysis (IPCA) on marker gene data and project cells onto principal components.

    Args:
        marker_mat (numpy.memmap):
            Memory-mapped array containing marker gene data (N-cells x len(marker_indices)).
        n_components_pca (int):
            Number of PCA components to project cells onto.
        sample_size (float):
            Fraction of N-cells to sample and run PCA on (0 < float <= 1).

    Returns:
        numpy.ndarray:
            Projected matrix containing cells projected onto principal components.

    Notes:
    ------
    This function applies Incremental Principal Component Analysis (IPCA) to marker gene data,
    scales the data using StandardScaler, and returns a matrix with cells projected onto
    the specified number of principal components.

    Example:
    --------
    >>> marker_mat = np.memmap("marker_matrix.npy", dtype='float32', mode='r', shape=(1000, 50))
    >>> n_components_pca = 10
    >>> sample_size = 0.5
    >>> projected_matrix = run_ipca(marker_mat, n_components_pca, sample_size)
    # Runs IPCA on marker gene data and projects cells onto principal components.
    """
    print('running pca...')
    scaler = StandardScaler()
    marker_mat = scaler.fit_transform(marker_mat)
    ipca = IncrementalPCA(n_components_pca, batch_size=int(sample_size*marker_mat.shape[0]))
    ipca.fit(marker_mat)
    projected_matrix = ipca.transform(marker_mat) 
    return projected_matrix

def run_gaussian_mixture(projected_matrix, sample_size, n_components_gaussian):
    """
    Estimate a Gaussian Mixture Model (GMM) based on data in PCA space and perform clustering.

    Args:
        projected_matrix (numpy.ndarray):
            2D array containing cell data projected onto principal components (N-cells x N-principal components).
        n_components_gaussian (int):
            Number of Gaussian components to model.
        sample_size (float):
            Fraction of N-cells to sample and use for estimating GMM parameters (0 < float <= 1).

    Returns:
        tuple:
            A tuple containing two elements:
            - predict (numpy.ndarray): Cluster labels assigned to cells.
            - scores (numpy.ndarray): Log likelihood scores for each cell.

    Notes:
    ------
    This function estimates a Gaussian Mixture Model (GMM) based on the data in PCA space, 
    samples a fraction of cells for parameter estimation, performs a grid search for GMM 
    hyperparameters, and returns cluster labels and log likelihood scores.

    Example:
    --------
    >>> projected_matrix = np.array([[1.2, 0.5], [0.8, 0.3], [2.0, 1.2]], dtype='float32')
    >>> sample_size = 0.5
    >>> n_components_gaussian = 3
    >>> predict, scores = run_gaussian_mixture(projected_matrix, sample_size, n_components_gaussian)
    # Estimates GMM based on PCA-space data and performs clustering.
    """
    print('estimating gaussian mixture model...')
    rand_index = np.random.choice(range(len(projected_matrix)), replace=False, size=int(sample_size*len(projected_matrix)))
    
    # first estimate N components and cov matrix
    X = projected_matrix[rand_index]
    
    ###############this code is from: https://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_selection.html#sphx-glr-auto-examples-mixture-plot-gmm-selection-py
    param_grid = {
        "n_components": range(11,16),
        "covariance_type": [ "full","spherical", "tied", "diag"]
    }
    print('gaussian gridsearch')
    
    grid_search = GridSearchCV(
        GaussianMixture(), param_grid=param_grid, scoring=gmm_bic_score, verbose=10
    )
    grid_search.fit(X)


    df = pd.DataFrame(grid_search.cv_results_)[
        ["param_n_components", "param_covariance_type", "mean_test_score"]
    ]
    df["mean_test_score"] = -df["mean_test_score"]
    df = df.rename(
        columns={
            "param_n_components": "Number of components",
            "param_covariance_type": "Type of covariance",
            "mean_test_score": "BIC score",
        }
    )
    ###############

    x = df.loc[np.argmin(df['BIC score'])]
    print(x)

    gm = GaussianMixture(n_components=x[0], covariance_type=x[1], random_state=0).fit(X)
    predict = gm.predict(projected_matrix)
    scores = gm.score_samples(projected_matrix)
    return predict, scores

def get_celltype_summary_score_per_cell(unique_celltype_name, marker_out, celltype_name):
    """
    Compute per-cell marker gene expression averages for each specified cell type label.

    Args:
        unique_celltype_name (numpy.ndarray):
            Enumeration of each cell type label being considered.
        marker_out (numpy.memmap):
            Memory-mapped array containing marker gene data (N-cells x len(marker_indices)).
        celltype_name (numpy.ndarray):
            Cell type label corresponding to each marker gene in marker_gene names (i.e., marker_out column names);
            length = len(marker_indices).

    Returns:
        numpy.ndarray:
            Array containing per-cell summary scores for each specified cell type label.

    Notes:
    ------
    This function computes per-cell marker gene expression averages for each specified cell type label.
    It uses marker gene data and cell type labels to calculate summary scores for each cell type.

    Example:
    --------
    >>> unique_celltype_name = np.array(['CellTypeA', 'CellTypeB'])
    >>> marker_out = np.memmap("marker_matrix.npy", dtype='float32', mode='r', shape=(1000, 50))
    >>> celltype_name = np.array(['CellTypeA', 'CellTypeB', 'CellTypeA', ...], dtype='<U10')
    >>> celltype_summary_score = get_celltype_summary_score_per_cell(unique_celltype_name, marker_out, celltype_name)
    # Computes per-cell summary scores for specified cell type labels.
    """
    print('computing celltype summary scores...')
    celltype_summary_score = np.empty((len(unique_celltype_name), len(marker_out)))
    for i in range(len(unique_celltype_name)):
        index = celltype_name==unique_celltype_name[i]
        celltype_summary_score[i] = np.mean(marker_out[:, index], axis=1)
    return celltype_summary_score
      
def get_celltype_annotations_per_group(unique_celltype_name, celltype_summary_score, predictions):
    """
    Assign cell types to each group based on marker genes and Gaussian mixture model predictions.

    Args:
        unique_celltype_name (numpy.ndarray):
            Enumeration of each cell type label being considered.
        celltype_summary_score (numpy.ndarray):
            Per-cell marker gene expression averages for each specified cell type label (N-labels x N-cells).
        predictions (numpy.ndarray):
            Assigned Gaussian component based on Gaussian mixture model (length = N-cells).

    Returns:
        tuple:
            A tuple containing two elements:
            - annotations_dic (dict): Dictionary mapping group identifiers to assigned cell type labels.
            - group_annotation_matrix (numpy.ndarray): Annotation matrix (N-groups x N-labels).

    Notes:
    ------
    This function assigns cell types to each group based on marker gene expression averages
    and Gaussian mixture model predictions. It returns a dictionary of group annotations
    and an annotation matrix.

    Example:
    --------
    >>> unique_celltype_name = np.array(['CellTypeA', 'CellTypeB'])
    >>> celltype_summary_score = np.array([[0.8, 0.2], [0.5, 0.6]], dtype='float32')
    >>> predictions = np.array([0, 1, 0, 1], dtype='int32')
    >>> annotations_dic, group_annotation_matrix = get_celltype_annotations_per_group(unique_celltype_name, celltype_summary_score, predictions)
    # Assigns cell types to groups based on marker genes and GMM predictions.
    """
    print('annotating cell types...')
    cellgroup = np.unique(predictions)
    group_annotation_matrix = np.empty((len(cellgroup), len(celltype_summary_score)))
    for i in range(len(cellgroup)):
        index = np.array(predictions)==cellgroup[i]
        group_annotation_matrix[i] = np.mean(celltype_summary_score[:,index], axis = 1)
    group_annotation_matrix = group_annotation_matrix/(np.sum(group_annotation_matrix, axis = 0)) # normalization
    index = np.argmax(group_annotation_matrix, axis = 1)
    annotations_dic = dict(zip(cellgroup, unique_celltype_name[index]))
    return annotations_dic, group_annotation_matrix
    
def assign_major_celltypes(marker_genes, celltype_dic, marker_out, predictions):
    """
    Assign major cell type labels to each cell based on marker genes and Gaussian mixture model predictions.

    Args:
        marker_genes (numpy.ndarray):
            Marker gene names (marker_out column names).
        celltype_dic (dict):
            Dictionary assigning each predicted Gaussian label to a cell type annotation.
        marker_out (numpy.memmap):
            Memory-mapped array containing marker gene data (N-cells x len(marker_indices)).
        predictions (numpy.ndarray):
            Assigned Gaussian component based on Gaussian mixture model (length = N-cells).

    Returns:
        numpy.ndarray:
            Array containing major cell type labels assigned to each cell.

    Notes:
    ------
    This function assigns major cell type labels to each cell based on marker genes, a cell type dictionary,
    and Gaussian mixture model predictions.

    Example:
    --------
    >>> marker_genes = np.array(['GeneA', 'GeneB', 'GeneC'], dtype='<U17')
    >>> celltype_dic = {'GeneA': 'CellTypeA', 'GeneB': 'CellTypeB'}
    >>> marker_out = np.memmap("marker_matrix.npy", dtype='float32', mode='r', shape=(1000, 50))
    >>> predictions = np.array([0, 1, 0, 1], dtype='int32')
    >>> annotations = assign_major_celltypes(marker_genes, celltype_dic, marker_out, predictions)
    # Assigns major cell types to each cell based on marker genes and GMM predictions.
    """
    celltype_name = np.array([celltype_dic[marker_genes[x]] for x in range(len(marker_genes))])
    unique_celltype_name = np.unique(celltype_name)
    celltype_summary_score = get_celltype_summary_score_per_cell(unique_celltype_name, marker_out, celltype_name)
    annotations_dic, summary_score = get_celltype_annotations_per_group(unique_celltype_name, celltype_summary_score, predictions)
    annotations = np.array([annotations_dic[predictions[i]] for i in range(len(predictions))])
    return annotations

def get_predictions(filtered_counts, marker_indices, marker_out, total_counts, n_components_pca, sample_size, n_components_gaussian):
    """
    Perform cell labeling based on a Gaussian Mixture Model (GMM) and principal component analysis (PCA).

    Args:
        filtered_counts (numpy.ndarray):
            2D array of counts post-cell filtering (N-cells x N-features).
        marker_indices (list):
            List of column indices indicating marker genes in the filtered_counts matrix.
        marker_out (numpy.memmap):
            Memory-mapped array for storing the normalized marker matrix (N-cells x len(marker_indices)).
        total_counts (numpy.ndarray):
            1D array of per-cell total counts, length = N-cells.
        n_components_pca (int):
            Number of PCA components to project cells onto.
        sample_size (float):
            Fraction of N-cells to sample for labeling (0 < float <= 1).
        n_components_gaussian (int):
            Number of Gaussian components to model.

    Returns:
        tuple:
            A tuple containing three elements:
            - projected_matrix (numpy.ndarray): Matrix with cells projected onto principal components.
            - predict (numpy.ndarray): Cluster labels assigned to cells.
            - scores (numpy.ndarray): Log likelihood scores for each cell.

    Notes:
    ------
    This function performs cell labeling based on a Gaussian Mixture Model (GMM) and principal component analysis (PCA).
    It computes the normalized marker matrix, projects cells onto PCs, and assigns cluster labels using GMM.

    Example:
    --------
    >>> filtered_counts = np.array([[10, 5, 3], [8, 4, 2]], dtype='int32')
    >>> marker_indices = [0, 2]
    >>> marker_out = np.memmap("marker_matrix.npy", dtype='float32', mode='w+', shape=(2, 2))
    >>> total_counts = np.array([18, 14], dtype='int32')
    >>> n_components_pca = 10
    >>> sample_size = 0.5
    >>> n_components_gaussian = 3
    >>> projected_matrix, predict, scores = get_predictions(filtered_counts, marker_indices, marker_out, total_counts, n_components_pca, sample_size, n_components_gaussian)
    # Performs cell labeling based on GMM and PCA.
    """
    print('getting normalized marker matrix...')
    with ProgressBar(total=filtered_counts.shape[0]) as numba_progress:
        get_marker_matrix_norm(filtered_counts, np.array(marker_indices), marker_out, total_counts, numba_progress)
    projected_matrix = run_ipca(marker_out, n_components_pca, sample_size)
    predict, scores = run_gaussian_mixture(projected_matrix, sample_size, n_components_gaussian)
    return projected_matrix, predict, scores
    
def get_major_annotations(marker_path, gene_names, filtered_counts, total_counts, sample_size, n_components_pca, n_components_gaussian, out_dir, infer_N_markers):
    """
    Return per-cell predicted cell type annotations based on marker genes, GMM, and PCA.

    Args:
        marker_path (str):
            Path to a CSV file containing marker gene annotations.
        gene_names (numpy.ndarray):
            1D array containing gene names with '<U17' data type (array length = N-features).
        filtered_counts (numpy.ndarray):
            2D array of counts post-cell filtering (N-cells x N-features).
        total_counts (numpy.ndarray):
            1D array of per-cell total counts, length = N-cells.
        sample_size (float):
            Fraction of N-cells to sample for labeling (0 < float <= 1).
        n_components_pca (int):
            Number of PCA components to project cells onto.
        n_components_gaussian (int):
            Number of Gaussian components to model.
        out_dir (str):
            Output directory where results will be saved.
        infer_N_markers (int or None):
            Number of markers to infer for reduction and annotation. If None, prior markers will be used.

    Returns:
        tuple:
            A tuple containing seven elements:
            - annotations (numpy.ndarray): Per-cell predicted cell type annotations.
            - marker_out (numpy.memmap): Memory-mapped array containing the normalized marker matrix (N-cells x N-markers).
            - projected_matrix (numpy.ndarray): Matrix with cells projected onto principal components.
            - predict (numpy.ndarray): Cluster labels assigned to cells.
            - scores (numpy.ndarray): Log likelihood scores for each cell.
            - marker_genes (numpy.ndarray): Marker gene names.
            - marker_indices (list): Column indices indicating marker genes in filtered_counts matrix.

    Notes:
    ------
    This function returns per-cell predicted cell type annotations based on marker genes, a Gaussian Mixture Model (GMM),
    and Principal Component Analysis (PCA). It can either use inferred markers for reduction and annotation or use prior markers.

    Example:
    --------
    >>> marker_path = "marker_annotations.csv"
    >>> gene_names = np.array(['GeneA', 'GeneB', 'GeneC'], dtype='<U17')
    >>> filtered_counts = np.array([[10, 5, 3], [8, 4, 2]], dtype='int32')
    >>> total_counts = np.array([18, 14], dtype='int32')
    >>> sample_size = 0.5
    >>> n_components_pca = 10
    >>> n_components_gaussian = 3
    >>> out_dir = "results"
    >>> infer_N_markers = 10
    >>> annotations, marker_out, projected_matrix, predict, scores, marker_genes, marker_indices = get_major_annotations(marker_path, gene_names, filtered_counts, total_counts, sample_size, n_components_pca, n_components_gaussian, out_dir, infer_N_markers)
    # Returns cell type annotations, marker matrix, and other results.
    """
    if infer_N_markers is not None:
        print('inferring markers for reduction & annotation...')
        marker_indices = infer_markers_for_reduction(infer_N_markers, filtered_counts)
        marker_out = np.memmap(out_dir + '/marker_matrix_norm.npy', mode='w+', shape=(filtered_counts.shape[0], len(marker_indices)), dtype='float')
        projected_matrix, predict, scores = get_predictions(filtered_counts, marker_indices, marker_out, total_counts, n_components_pca, sample_size, n_components_gaussian)
        marker_genes = gene_names[marker_indices]
        out_markers, out_fc = infer_marker_genes_per_cluster(marker_out, predict, marker_genes)
        annotations = np.concatenate((out_markers, out_fc))
    else:
        print('using prior markers for reduction & annotation...')
        marker_indices, celltype_dic = get_marker_indices(marker_path, gene_names)
        marker_out = np.memmap(out_dir + '/marker_matrix_norm.npy', mode='w+', shape=(filtered_counts.shape[0], len(marker_indices)), dtype='float')
        projected_matrix, predict, scores = get_predictions(filtered_counts, marker_indices, marker_out, total_counts, n_components_pca, sample_size, n_components_gaussian)
        marker_genes = gene_names[marker_indices]
        annotations = assign_major_celltypes(marker_genes, celltype_dic, marker_out, predict)
    return annotations, marker_out, projected_matrix, predict, scores, marker_genes, marker_indices

def infer_markers_for_reduction(infer_N_markers, filtered_counts):
    """
    Infer marker genes for dimensionality reduction and annotation based on gene expression variance.

    Args:
        infer_N_markers (int):
            Number of markers to infer for dimensionality reduction and annotation.
        filtered_counts (numpy.ndarray):
            2D array of counts post-cell filtering (N-cells x N-features).

    Returns:
        numpy.ndarray:
            Array containing the indices of inferred marker genes.

    Notes:
    ------
    This function infers marker genes for dimensionality reduction and annotation based on gene expression trends.
    It identifies genes with high mean-variance trends and selects the top N markers for further analysis.

    Example:
    --------
    >>> infer_N_markers = 10
    >>> filtered_counts = np.array([[10, 5, 3], [8, 4, 2]], dtype='int32')
    >>> marker_indices = infer_markers_for_reduction(infer_N_markers, filtered_counts)
    # Infers marker genes for dimensionality reduction and annotation.
    """
    mean = get_log2(np.mean(filtered_counts, axis = 0))
    var = get_log2(np.var(filtered_counts, axis = 0))
    var_predictions = localreg.localreg(mean, var)
    mean_var_trend_removed = var-var_predictions
    marker_indices = np.argsort(-1*(mean_var_trend_removed + mean))[0:infer_N_markers]
    return marker_indices

def infer_marker_genes_per_cluster(marker_out, predict, marker_genes):
    """
    Infer marker genes per cluster based on fold-change (FC).

    Args:
        marker_out (numpy.memmap):
            Memory-mapped array containing the normalized marker matrix (N-cells x N-markers).
        predict (numpy.ndarray):
            Cluster labels assigned to cells.
        marker_genes (numpy.ndarray):
            Marker gene names.

    Returns:
        tuple:
            A tuple containing two elements:
            - out_markers (numpy.ndarray): Marker gene names per cluster (N-clusters x n_markers).
            - out_fc (numpy.ndarray): Fold-change values per cluster (N-clusters x n_markers).

    Notes:
    ------
    This function infers marker genes per cluster based on fold-change (FC).
    For each cluster, it computes the FC of marker genes for in-cluster cells compared to out-cluster cells and selects the top N markers.

    Example:
    --------
    >>> marker_out = np.memmap("marker_matrix.npy", dtype='float32', mode='r', shape=(1000, 50))
    >>> predict = np.array([0, 1, 0, 1], dtype='int32')
    >>> marker_genes = np.array(['GeneA', 'GeneB', 'GeneC'], dtype='<U17')
    >>> out_markers, out_fc = infer_marker_genes_per_cluster(marker_out, predict, marker_genes)
    # Infers marker genes per cluster based on FC analysis.
    """    
    n_markers = 20
    unique_predictions = np.unique(predict)
    out_fc = np.empty(shape=(len(unique_predictions), n_markers), dtype='float')
    out_markers = np.empty(shape=(len(unique_predictions), n_markers), dtype='<U17')
    for i in range(len(unique_predictions)):
        current_prediction = unique_predictions[i]
        FC = np.log2(np.mean(marker_out[np.where(predict==current_prediction)[0]], axis = 0)/np.mean(marker_out[np.where(predict!=current_prediction)[0]], axis = 0))
        indices = np.argsort(-1*FC)[0:n_markers]
        out_markers[i] = marker_genes[indices]
        out_fc[i] = FC[indices]
    return out_markers, out_fc

def get_log2(vals):
    """
    Compute the logarithm base 2 of non-zero values after shifting.

    Args:
        vals (numpy.ndarray):
            1D array containing numerical values.

    Returns:
        numpy.ndarray:
            1D array containing the logarithm base 2 of non-zero values after shifting.

    Notes:
    ------
    This function computes the logarithm base 2 of non-zero values in the input array after shifting the values.
    Shifting ensures that zero values are excluded from the logarithmic operation.

    Example:
    --------
    >>> vals = np.array([0, 1, 2, 3], dtype='float32')
    >>> log_vals = get_log2(vals)
    # Computes the logarithm base 2 of non-zero values after shifting.
    """  
    vals+=np.min(vals[vals!=0])
    vals = np.log2(vals)
    return vals
