import matplotlib.backends.backend_pdf
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE

def plot_filter_cells(mito_fractions, total_counts, keep_cells, pdf, total_counts_lower_bound, total_counts_upper_bound):
    data = pd.DataFrame(np.hstack((mito_fractions.reshape(-1,1), (total_counts.reshape(-1,1)), keep_cells.reshape(-1,1))))
    data.columns = ['mito_fractions', 'total_counts', 'group']
    #data = data[(data['total_counts']>total_counts_lower_bound) & (data['total_counts']<total_counts_upper_bound)]
    sns.jointplot(x="mito_fractions", y="total_counts", data=data, hue='group', joint_kws=dict(alpha=0.5), s = 2)
    pdf.savefig()
    pdf.close()

def plot_cells_reduced_dim(projected_matrix, annotations):
    X_embedded = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(projected_matrix)
    data = pd.DataFrame(np.hstack((annotations.reshape(-1,1), X_embedded[:,0:2])))
    data.columns = ['anno', 'DIM1', 'DIM2']
    data['DIM1'] = pd.to_numeric(data['DIM1'])
    data['DIM2'] = pd.to_numeric(data['DIM2'])
    sns.jointplot(x="DIM1", y="DIM2", data=data, hue='anno', joint_kws=dict(alpha=0.5))