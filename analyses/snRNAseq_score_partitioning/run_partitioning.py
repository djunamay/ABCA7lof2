# run many KL and METIS iterations on the Ex-LE graph

import numpy as np
import gseapy
from scipy.sparse import csr_matrix
from tqdm.contrib.concurrent import process_map
from functools import partial
from scipy.sparse.csgraph import shortest_path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import networkx
import numba as nb
from sklearn.manifold import SpectralEmbedding
from tqdm import tqdm
from sklearn.cluster import SpectralClustering
import seaborn as sns
import metis
import networkx as nx
import metis

from ABCA7lof2.geneclusters import evaluate_cut, get_scores, get_kernighan_lin_clusters, get_gene_pathway_matrix, get_full_matrix_from_bipartite, plot_component, plot_edges, plot_nodes, group, compute_groupped_matrix, get_scores, find_similar_clusters, get_representative_name_per_cluster, get_kernighan_lin_clusters, get_gene_pathway_matrix, compute_groupped_matrix, get_full_matrix_from_bipartite

import os
import urllib.request

output_dir = "../../processed_data/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

urls = [
    'https://storage.googleapis.com/abca7lof/scRNAseq/processed_data/leading_edge_0825Ex.csv',
    'https://storage.googleapis.com/abca7lof/scRNAseq/processed_data/WikiPathways_2019_Human.npy'
]
for url in urls:
    local_filename = os.path.join(output_dir, os.path.basename(url))
    urllib.request.urlretrieve(url, local_filename)
    print("Download completed!")

# save wikipaths graph
p1 = np.load('./raw_data/genesets/WikiPathways_2019_Human.npy', allow_pickle=True).item()
res = {**p1}
np.save('./processed_data/genesets/all_paths.npy', res)
pd.DataFrame.from_dict(res, orient='index').to_csv('./processed_data/genesets/all_paths.csv')

# subset matrix
mat = get_gene_pathway_matrix('./processed_data/genesets/all_paths.npy')
leading_edge = './processed_data/for_plotting/leading_edge_0825Ex.csv'
leading_edge = pd.read_csv(leading_edge, index_col=0)
S = set(leading_edge['gene'])

col_index = np.where([x in S for x in mat.columns])[0]
mat_sub = mat.iloc[:,col_index]

path_index = (np.sum(mat_sub, axis=1)>4)
mat_sub = mat_sub.loc[path_index]
mat_sub = mat_sub.loc[:,np.sum(mat_sub, axis=0)>0]

# run many KL iterations
C = 0
KL_modified = True
random_labels = True
unweighted = True

N=50000
loss_kl = np.empty(N)
labs_kl = np.empty((N,np.sum(mat_sub.shape)))

for i in tqdm(range(N)):
    frame, loss_temp = get_kernighan_lin_clusters(None, 50, C, KL_modified, random_labels, unweighted, seed=i, no_progress=True, mat=mat_sub)
    frame.columns = ['cluster', 'description', 'is_gene']
    labs_kl[i] = np.array(frame['cluster'])
    loss_kl[i] = loss_temp
    
np.save('./processed_data/genesets/kl_labs.npy', labs_kl)
np.save('./processed_data/genesets/kl_loss.npy', loss_kl)

# compute laplacian matrix
full_mat = np.zeros((np.sum(mat_sub.shape), np.sum(mat_sub.shape)))
full_mat[mat_sub.shape[1]:][:,:mat_sub.shape[1]] = mat_sub
full_mat[:mat_sub.shape[1]][:,mat_sub.shape[1]:]=mat_sub.T

g = nx.from_numpy_matrix(full_mat)
clusters = np.empty((N,np.sum(mat_sub.shape)))

# run METIS iterations
loss_met=np.empty(N)
labs_met = np.zeros_like(clusters)
for i in tqdm(range(N)): 
    sc = metis.part_graph(g, nparts=8, tpwgts=None, ubvec=None, recursive=False, seed=i)[1]
    loss_met[i]=evaluate_cut(np.ascontiguousarray(mat_sub.values.T), sc, 0)
    labs_met[i] = sc

np.save('./processed_data/genesets/met_labs.npy', labs_met)
np.save('./processed_data/genesets/met_loss.npy', loss_met)