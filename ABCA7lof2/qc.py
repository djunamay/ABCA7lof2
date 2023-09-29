import numba as nb
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.mixture import GaussianMixture
from numba_progress import ProgressBar
import ipdb 
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GridSearchCV
import pandas as pd

@nb.njit(parallel=True)
def filter_logcounts(temp_output, logcounts, index_cells, index_genes):
    for i in nb.prange(len(index_cells)):
        temp_output[i] = logcounts[index_cells[i]][index_genes]

def compute_thresh(adata, obs_slot, grp_slot):
    df = adata.obs[[obs_slot, grp_slot]].groupby(grp_slot).mean()
    x = df[obs_slot]
    thresh_up = np.mean(x) + 1.5*np.std(x)
    thresh_down = np.mean(x) - 1.5*np.std(x)
    return df, thresh_up, thresh_down

def filter_clusters(temp, grp, remove):
    df, up, down = compute_thresh(temp, 'mito_fractions', grp)
    x = list(np.array(df.index)[np.array(df['mito_fractions'])>up])
    for i in x:
        remove.append(i)

    df, up, down = compute_thresh(temp, 'total_counts', grp)
    x = list(np.array(df.index)[(np.array(df['total_counts'])>up)]) 
    for i in x:
        remove.append(i)
    
    i = np.argsort(np.array(temp.obs[grp]))
    temp = temp[i]

    cols = [x in set(remove)for x in np.array(temp.obs[grp])]

    sns.set(rc={'figure.figsize':(20,3)})

    sns.boxplot(x=np.array(temp.obs[grp]), y=np.array(temp.obs['mito_fractions']), hue=cols)
    plt.ylabel('mito_fractions')
    plt.legend([], [], frameon=False)

    plt.show()

    sns.boxplot(x=np.array(temp.obs[grp]), y=np.array(temp.obs['total_counts']), hue=cols)
    plt.ylabel('total_counts')
    plt.legend([], [], frameon=False)

    plt.show()

    b = pd.DataFrame(np.unique(temp.obs[grp], return_counts=True)).T
    b.columns = ['cluster', 'N']
    b['COL'] =[ x in set(remove) for x in b['cluster'] ]
    sns.barplot(data = b, x = 'cluster', y = 'N', hue = 'COL')
    plt.ylabel('N cells')
    plt.legend([], [], frameon=False)

    plt.show()

    unique_samples = list()
    temp2 = temp.obs[[grp, 'sample_id']]
    for i in np.unique(temp2[grp]):
        unique_samples.append(len(np.unique(temp2['sample_id'][temp2[grp]==i])))
    COL = [ x in set(remove) for x in np.unique(temp2[grp])]
    sns.barplot(x = np.unique(temp2[grp]), y = unique_samples, hue = COL)
    plt.ylabel('N individuals')
    plt.legend([], [], frameon=False)


    sns.set(rc={'figure.figsize':(5,3.5)})
    temp.obs['COL'] = np.array([x if x in set(remove) else 'other' for x in temp.obs[grp]])
    sc.pl.umap(temp, color=['COL'], add_outline=False, frameon=False, gene_symbols='Gene', return_fig=False)

@nb.njit(parallel=False)
def filter_genes(counts, keep_genes):
    i = counts.shape[0]
    for x in range(i):
        keep_genes |= counts[x] > 0

@nb.njit(parallel=True)
def get_total_counts(counts, total_counts):
    for i in nb.prange(counts.shape[0]):
        total_counts[i] = counts[i].sum()
    return total_counts

@nb.njit(parallel=True)
def compute_mito_fraction(counts, gene_index, mito_fractions, total_counts):
    for i in nb.prange(counts.shape[0]):
        total_counts[i] = counts[i].sum()
        if total_counts[i]>0:
            mito_fractions[i] = counts[i][gene_index].sum() / total_counts[i]  # .sum()
        else:
            mito_fractions[i] = counts[i][gene_index].sum() / (total_counts[i]+1)

def compute_total_counts(counts, total_counts):
    for i in range(counts.shape[0]):
        total_counts[i] = counts[i].sum()
    
def get_fraction_mito(gene_names, counts, prefix):
    '''
    returns per-cell mitochondrial fraction estimates
    returns per-cell total counts (sum)
    returns gene-wise index of mitochondrial genes
    Args:
        gene_names: ndarray
            1D array containing gene names with '<U17' data type; array length = N-features
        prefix: string
            mitochondrial genes start with this prefix (e.g. 'MT-')
        counts: ndarray
            2D array containing data with 'int16' data type
            N-cells x N-features
    '''
    mito_index = np.char.startswith(gene_names, prefix)
    mito_fractions = np.empty(counts.shape[0])
    total_counts = np.empty(counts.shape[0], dtype=int)
    compute_mito_fraction(counts, mito_index, mito_fractions, total_counts)
    return mito_fractions, mito_index, total_counts


def cluster_on_mito_and_counts(mito_fractions, total_counts, sample_size):
    '''
    returns boolean array indicating which cells to keep based on gaussian mixture model on total counts and fraction mitochondrial reads
    Args:
        mito_fractions: ndarray
            1D array containing per-cell mitochondrial fractions (dtype = 'float')
        total_counts: ndarray
            1D array containing per-cell total counts (dtype = 'int')
        sample_size: float
            0< float <=1 indicating the fraction of N-cells to sample uniformly to estimate parameters of gaussian mixture model
    '''
    N = len(mito_fractions)
    rand_index = np.random.choice(range(N), replace=False, size=int(sample_size * N))
    complete_data = np.stack([mito_fractions, total_counts], axis=1)
    sub_sampled = complete_data[rand_index]
    #ipdb.set_trace()
    gm = GaussianMixture(n_components=2, random_state=0).fit(sub_sampled)
    labels = gm.predict(complete_data)
    keep_cells = labels == np.bincount(labels).argmax()
    stats = np.unique(labels, return_counts=True)
    print('Labels = ' + str(stats[0]))
    print('counts per label = ' + str(stats[1]))
    return keep_cells

def filter_cells_by_mito(mito_fractions):
    
    X  = np.log10(mito_fractions+np.min(mito_fractions[mito_fractions>0])).reshape(-1,1)
    ###############this code is from: https://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_selection.html#sphx-glr-auto-examples-mixture-plot-gmm-selection-py
    param_grid = {
        "n_components": range(1, 15),
        "covariance_type": ["spherical", "tied", "diag", "full"],
    }
    print('gaussian gridsearch')
    
    grid_search = GridSearchCV(
        GaussianMixture(), param_grid=param_grid, scoring=gmm_bic_score
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
    
    gm = GaussianMixture(n_components=x[0], max_iter=10000, n_init=100, covariance_type=x[1]).fit(X)
    labels = gm.predict(X)
    clust_remove = np.argmax([np.mean(X[labels==x]) for x in np.unique(labels)])

    keep_cells = labels!=clust_remove
    return keep_cells

def filter_cells(gene_names, counts, meta, sample_size, total_counts_lower_bound, total_counts_upper_bound, pdf, prefix, out_dir):
    '''
    filters counts matrix based on per-cell total counts and per-cell mitochondrial fraction
    Args:
        gene_names ndarray
            1D array containing gene names with '<U17' data type; array length = N-features
        counts ndarray
            2D array containing data with 'int16' data type
            N-cells x N-features
        meta ndarray
            2D array, N-cells x N-metadata-features
        sample_size: float
            0< float <=1 indicating the fraction of N-cells to sample uniformly to estimate parameters of gaussian mixture model
        total_counts_lower_bound integer
            minimum number of total counts per cell for a cell to be considered
        total_counts_upper_bound
            maximum number of total counts per cell for a cell to be considered
        pdf A PDF Matplotlib backend
        prefix: string
            mitochondrial genes start with this prefix (e.g. 'MT-')
    '''
    print('getting mito fractions')
    mito_fractions, mito_index, total_counts = get_fraction_mito(gene_names, counts, prefix)  
    
    print('searching for best params')
    keep_cells = filter_cells_by_mito(mito_fractions)

    print('All cells')
    print(len(keep_cells))
    print('Mito keep')
    print(np.sum(keep_cells)/len(keep_cells))
    
    with ProgressBar(total=counts.shape[0]) as numba_progress:
        N = np.empty(len(counts))
        get_N_detected_genes(counts, N, numba_progress)
    
    keep_cells = keep_cells&(N>=total_counts_lower_bound)&(N<=total_counts_upper_bound)
    
    print('All keep:')
    print(np.sum(keep_cells)/len(keep_cells))
    
    filtered_counts = np.memmap(out_dir + '/filtered_counts.npy', mode='w+', shape=(np.count_nonzero(keep_cells), counts.shape[1]), dtype='int16')
    meta_subset = meta[keep_cells]
    mito_fractions_subset = mito_fractions[keep_cells]

    print('Mito fraction mean discard:')
    print(np.mean(mito_fractions[np.invert(keep_cells)]))
    
    print('Mito fraction mean keep:')
    print(np.mean(mito_fractions_subset))
    
    total_counts_subset = total_counts[keep_cells]
    with ProgressBar(total=filtered_counts.shape[0]) as numba_progress:
        filter_counts_by_keep_cells(np.where(keep_cells)[0], counts, filtered_counts, numba_progress)
    return mito_index, total_counts_subset, mito_fractions_subset, meta_subset, filtered_counts

############### this function is from: https://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_selection.html#sphx-glr-auto-examples-mixture-plot-gmm-selection-py
def gmm_bic_score(estimator, X):
    """Callable to pass to GridSearchCV that will use the BIC score."""
    # Make it negative since GridSearchCV expects a score to maximize
    return -estimator.bic(X)
###############


def filter_cells_by_major_annotation(mito_fractions, total_counts, sample_size, celltype_annotations, individual_annotation, filtered_counts):
    print('Filtering Cells by major annotation.')
    keep_cells_0 = filter_cells_by_celltype(mito_fractions, total_counts, sample_size, celltype_annotations)
    keep_cells_1, keep_ind = filter_individuals(celltype_annotations, individual_annotation)
    keep_genes_per_celltype = filter_genes_by_expression(celltype_annotations, filtered_counts)
    return keep_cells_0 + keep_cells_1, keep_cells_0, keep_cells_1, keep_genes_per_celltype, keep_ind


def filter_cells_by_celltype(mito_fractions, total_counts, sample_size, annotation):
    unique_celltype = np.unique(annotation)
    keep_cells_all = np.empty(len(annotation)).astype('bool')
    for i in range(len(unique_celltype)):
        index = annotation==unique_celltype[i]
        keep_cells_all[index] = cluster_on_mito_and_counts(mito_fractions[index], total_counts[index], sample_size)
    return keep_cells_all.astype('bool') # will this end up removing too many cells? --> have to check

def filter_individuals(celltype_annotations, individual_annotation):
    table = pd.crosstab(celltype_annotations, individual_annotation)
    individual_names = table.columns
    table = np.array(table)
    fractions = table / np.sum(table, axis=0)
    means = np.mean(fractions, axis=1)
    stds = np.std(fractions, axis=1) * 2
    discard_individuals = np.abs(fractions - means.reshape(-1, 1)) > stds.reshape(-1, 1)
    keep_individuals = np.array(individual_names)[np.sum(discard_individuals, axis=0) == 0]
    keep_cells = np.isin(individual_annotation, keep_individuals)
    return keep_cells.astype('bool'), keep_individuals

def filter_genes_by_expression(annotations, counts):
    celltypes = np.unique(annotations)
    N = counts.shape[1]
    keep_genes_all = np.empty((len(celltypes), N), dtype=bool)
    for i in range(len(celltypes)):
        index = annotations == celltypes[i]
        with ProgressBar(total=N) as numba_progress:
            filter_genes_by_celltype(N, index, counts, keep_genes_all, i, numba_progress)
    return keep_genes_all#.astype('bool')

@nb.njit(parallel=True)
def filter_genes_by_celltype(N, index, counts, keep_genes_all, i, progress_hook):
    for x in nb.prange(N):
        #ipdb.set_trace()
        gene = counts[:,x][index]
        keep_genes_all[i, x] = (np.sum(gene>0)/len(gene)) > 0.1
        progress_hook.update(1)

@nb.njit(parallel=True)
def filter_counts_by_keep_cells(counts_index, counts, filtered_counts, progress_hook):
    for i in nb.prange(len(counts_index)):
        filtered_counts[i] = counts[counts_index[i]]
        progress_hook.update(1)

@nb.njit(parallel=True)
def get_N_detected_genes(counts, N, progress_hook):
    for i in nb.prange(len(counts)):
        N[i] = sum(counts[i]>0)
        progress_hook.update(1)
        
def filter_on_gaussian_logliklihood(scores):
    scores = scores.reshape(-1,1)
    gm = GaussianMixture(n_components=3, covariance_type='full', random_state=0).fit(scores)
    predict = gm.predict(scores)
    keep_cells = predict!=np.bincount(predict).argmin()
    return keep_cells
    
def log_normalize_counts(counts, total_counts, norm_counts):
    median_counts = np.median(total_counts)
    for i in tqdm(range(counts.shape[0])):
        norm = median_counts*(counts[i]/total_counts[i])
        norm_counts[i] = np.log1p(norm, where=norm!=0, out=np.zeros_like(norm).astype('float'))