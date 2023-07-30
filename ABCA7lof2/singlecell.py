import numba as nb
from numba.typed import Dict

@nb.njit(parallel=True)
def bulk_data(logcounts, exclude, dictionary1, dictionary2, bulked, sums, projids, celltype):
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
    for i in nb.prange(counts.shape[0]):
        if i in exclude:
            continue
        else:
            expressed[celltype[i]] += counts[i]>0
            summed[celltype[i]] += 1
    return summed, expressed