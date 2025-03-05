import sys
sys.path.append('/Users/djuna/Documents/ABCA7lof2/')

import numpy as np
from tqdm import tqdm
import json
from ABCA7lof2.geneclusters import get_gene_pathway_matrix, get_kernighan_lin_clusters

import argparse

# Initialize the parser
parser = argparse.ArgumentParser(description='Assign leading edge genes to clusters based on gene-pathway graph')
parser.add_argument('--name', type=str, required=True, help='Path to the input JSON file containing leading edge genes')
args = parser.parse_args()
name = args.name

# load the leading edge genes
with open("../../processed_data/bulkRNAseq_fgsea_leading_edge.json", 'r') as file:
    data = json.load(file)

mat = get_gene_pathway_matrix('../../processed_data/all_paths.npy')
S = set(data[name])

# prep data
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

np.save('../../processed_data/kl_labs_' + name + '.npy', labs_kl)
np.save('../../processed_data/kl_loss_' + name + '.npy', loss_kl)