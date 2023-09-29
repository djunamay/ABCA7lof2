import numba as nb
from numba.typed import Dict

@nb.njit(parallel=True)
def bulk_data(logcounts, exclude, dictionary1, dictionary2, bulked, sums, projids, celltype):
    '''
    Process logcounts and update bulked and sums arrays for each cell based on specified criteria.

    Args:
        logcounts (numpy.ndarray): A 2D array of log-transformed counts for individual cells.
        exclude (list): A list of indices to exclude from processing.
        dictionary1 (dict): A dictionary mapping projids to indices.
        dictionary2 (dict): A dictionary mapping cell types to indices.
        bulked (numpy.ndarray): A 1D array to store bulked data.
        sums (numpy.ndarray): A 1D array to store sums of processed data.
        projids (numpy.ndarray): An array containing projids for each cell.
        celltype (numpy.ndarray): An array containing cell types for each cell.

    Returns:
        bulked (numpy.ndarray): Updated bulked data.
        sums (numpy.ndarray): Updated sums.
    '''
    for i in nb.prange(logcounts.shape[0]):
        if i in exclude:
            continue
        else:
            index = dictionary1[projids[i]]+dictionary2[celltype[i]]
            bulked[index]+=logcounts[i]
            sums[index]+=1
    return bulked, sums

@nb.njit(parallel=True)
def expressed_fraction(counts, exclude, celltype, expressed, summed):
    '''
    Calculate the fraction of expressed genes per cell type based on counts data.

    Args:
        counts (numpy.ndarray): A 1D array of gene expression counts for individual cells.
        exclude (list): A list of indices to exclude from processing.
        celltype (numpy.ndarray): An array containing cell types for each cell.
        expressed (numpy.ndarray): An array to store the count of expressed genes per cell type.
        summed (numpy.ndarray): An array to store the total count of genes per cell type.

    Returns:
        summed (numpy.ndarray): Updated total count of genes per cell type.
        expressed (numpy.ndarray): Updated count of expressed genes per cell type.
    '''
    for i in nb.prange(counts.shape[0]):
        if i in exclude:
            continue
        else:
            expressed[celltype[i]] += counts[i]>0
            summed[celltype[i]] += 1
    return summed, expressed