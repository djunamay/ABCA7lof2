import sys
sys.path.append('/Users/djuna/Documents/ABCA7lof2/')

from tqdm import tqdm
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
import json

from ABCA7lof2.plotting_geneclusters import get_le_clusters, get_layout, plot_graph, plot_rep_names, get_top_genes
from ABCA7lof2.geneclusters import get_gene_pathway_matrix, compute_jaccard_all_clust
from ABCA7lof2.geneclusters import evaluate_cut, get_scores, get_kernighan_lin_clusters, get_gene_pathway_matrix, get_full_matrix_from_bipartite, plot_component, plot_edges, plot_nodes, group, compute_groupped_matrix, get_scores, find_similar_clusters, get_representative_name_per_cluster, get_kernighan_lin_clusters, get_gene_pathway_matrix, compute_groupped_matrix, get_full_matrix_from_bipartite

from scipy.optimize import linear_sum_assignment
from statsmodels.stats.multitest import multipletests
import matplotlib.patches as mpatches
import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# Add custom legend; blue square and a red square
import matplotlib.patches as mpatches

def get_representative_name_per_cluster(bipartite_mat, colnames_mat, rownames_mat, description_table, cluster, cluster_column, N=5):
    genes = set(description_table.loc[description_table['is_gene']&(description_table[cluster_column]==cluster)]['description'])
    paths = set(description_table.loc[np.invert(description_table['is_gene'])&(description_table[cluster_column]==cluster)]['description'])
    if len(paths)==0:
        return 'C.'+str(cluster), np.nan, np.nan, np.nan
    else:
        index_col = [x in genes for x in colnames_mat]
        index_row = [x in paths for x in rownames_mat]
        
        #ipdb.set_trace()
    
        sum_internal = np.sum(bipartite_mat[index_row][:,index_col], axis=1)
        sum_external = np.sum(bipartite_mat[index_row][:,np.invert(index_col)], axis=1)
        sum_ratio = sum_internal#/(sum_external+sum_internal)
        S = np.argsort(-1*sum_ratio)[:N]
        rep_name = rownames_mat[index_row][S]

        return rep_name #'C.'+str(cluster), rep_name.split(' (')[0], sum_internal[S], sum_internal[S]

# get rep names

def plot_rep_names(pos, unique_clusters, colors, mat_sub, frame, out_path, N, cur_labels):
    import re
    from adjustText import adjust_text
    colnames = np.array(mat_sub.columns)
    rownames = np.array(mat_sub.index)

    out = [get_representative_name_per_cluster(np.array(mat_sub), colnames, rownames, frame, x, N) for x in np.unique(frame['cluster'])]

    plt.figure(figsize = (1,1))
    a = plt.gca()
    a.axis('off')
    texts = []
    y = 0
    for i, cluster_name in enumerate(unique_clusters):
            index = cur_labels==cluster_name
            x, y = np.mean(pos[index], axis=0)

            props = dict(boxstyle='round', facecolor='white', alpha=1, edgecolor=colors[i])
            temp = [re.split(r' [(]GO|WP| R-',x)[0] for x in out[i]]
            T = ('\n').join(temp)
            plt.text(0, 0, T,  bbox=props, c ='black', fontsize = 10)#, #style = "italic")#colors[i]

            a = plt.gca()
            a.axis('off')
            plt.savefig(out_path + str(i) + '.pdf', bbox_inches="tight")
            plt.figure()
            #y+=0.01

def return_partitioning(pathways_path, celltype, kl_loss_path, data):
    mat = get_gene_pathway_matrix(pathways_path)
    kl_loss = np.load(kl_loss_path)
    seed = np.argmin(kl_loss)
    genes = set(data[celltype])
    frame, mat_sub = get_le_clusters(None, None, mat, seed, 50, S=genes)

    return frame, mat_sub

def get_random_overlap(frame1, frame2, seed, size1, size2):
    x = set(frame1['description'].sample(n=size1, replace=False, random_state=seed))
    y = set(frame2['description'].sample(n=size2, replace=False, random_state=seed))
    return len(x.intersection(y))/len(x.union(y))

def get_true_overlap(frame1, frame2):
    x = set(frame1['description'])
    y = set(frame2['description'])
    return len(x.intersection(y))/len(x.union(y))

def adjust_pvalues(pvalues):
    pvalues_flat = pvalues.flatten()
    _, pvalues_corrected_flat, _, _ = multipletests(pvalues_flat, method='fdr_bh')
    pvalues_corrected = pvalues_corrected_flat.reshape(pvalues.shape)
    return pvalues_corrected

def sort_matrix(M):
    cost_matrix = -M
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    sorted_assignment = sorted(zip(col_ind, row_ind))  # sort by column index
    sorted_cols, sorted_rows = zip(*sorted_assignment)
    assigned_rows = list(sorted_rows)
    remaining_rows = [r for r in range(M.shape[0]) if r not in assigned_rows]
    new_row_order = assigned_rows + remaining_rows
    new_col_order = list(range(M.shape[1]))
    return new_row_order, new_col_order

def get_jaccard_matrix(frame1, frame2, iterations=10000):

    jaccard = np.zeros((len(np.unique(frame1['cluster'])), len(np.unique(frame2['cluster']))))
    pvalues = jaccard.copy()

    for i in tqdm(range(jaccard.shape[0])):
        frame1_temp = frame1[frame1['cluster']==i]
        for j in range(jaccard.shape[1]):
            frame2_temp = frame2[frame2['cluster']==j]

            true_overlap = get_true_overlap(frame1_temp, frame2_temp)
            random_overlaps = [get_random_overlap(frame1, frame2, seed, len(frame1_temp), len(frame2_temp)) for seed in range(iterations)]
            p_value = np.sum(np.array(random_overlaps) >= true_overlap)/iterations

            jaccard[i,j] = true_overlap
            pvalues[i,j] = p_value

    return jaccard, pvalues

def plot_jaccard_matrix(jaccard, adjusted_pvalues, new_row_order, new_col_order, title, xlabel, ylabel, label_y=-0.5):
    nrows, ncols = jaccard.shape
    jaccard_matrix = jaccard[new_row_order, :][:, new_col_order]
    pvalue_matrix = adjusted_pvalues[new_row_order, :][:, new_col_order]

    im = plt.imshow(jaccard_matrix, cmap='viridis', )
    plt.title(title, fontsize=10, pad=10)
    # Annotate each cell with its p-value
    for i in range(nrows):
        for j in range(ncols):
            # Format the p-value with 3 decimal places
            if pvalue_matrix[i, j] < 0.05:
                p_text = "*"
                plt.text(j, i, p_text, ha="center", va="center", color="w", fontsize=15)

    plt.xlabel(xlabel, labelpad=10)
    plt.ylabel(ylabel, labelpad=10)

    plt.legend().set_visible(False)
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.set_label('Jaccard Index')
    # Custom legend

    # Create a custom legend
    label = r'* = empirical $p_{adj} < 0.05$'
    #label = r'\textcolor{red}{*} = empirical $p_{adj} < 0.05$'

    #label = r'$\color{red}{*}$ = empirical $p_{adj} < 0.05$'

    legend_elements = [mpatches.Patch(color='white', label=label)]#, mpatches.Patch(color='red', label='up-regulated cluster'), mpatches.Patch(color='blue', label='down-regulated cluster')]
    plt.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.25, label_y), frameon=False)


def get_scores_per_cluster(data, frame, cluster_column):
    g_names  = np.array(frame['description'][frame['is_gene']])
    scores = pd.DataFrame(data['score'])
    scores.columns = ['iN']
    SCORES = np.array(scores.loc[g_names])

    names = np.unique(frame[cluster_column])
    temp = [SCORES[frame[cluster_column][frame['is_gene']]==i] for i in names]

    return dict(zip(names, temp))

def get_names_per_cluster(mat_sub, frame, cluster_column, N=5):
    out = [get_representative_name_per_cluster(np.array(mat_sub), np.array(mat_sub.columns), np.array(mat_sub.index), frame, x, cluster_column, N) for x in np.unique(frame[cluster_column])]
    return out

def plot_direction_colors(cluster_labels_1, cluster_labels_2, scores_1, scores_2, new_row_order, new_col_order, ax):

    g2_labels = np.array(cluster_labels_1)[new_row_order]
    y_labels = np.array(cluster_labels_2)[new_col_order]
    g2_means = ['down' if np.mean(scores_1[int(ix.split('.')[1])])<0 else 'up' for ix in g2_labels]
    y_means = ['down' if np.mean(scores_2[int(ix.split('.')[1])])<0 else 'up' for ix in y_labels]

    for tick_label, direction in zip(ax.get_xticklabels(), y_means):
        tick_label.set_color('red' if direction == 'up' else 'blue')

    # Color x tick labels similarly
    for tick_label, direction in zip(ax.get_yticklabels(), g2_means):
        tick_label.set_color('red' if direction == 'up' else 'blue')

def map_clusters_rows_to_cols(jaccard, adjusted_pvalues, row_prefix):

    binary_pvalues = adjusted_pvalues<0.05
    cluster_assignment = np.zeros(adjusted_pvalues.shape[0])

    for i in range(adjusted_pvalues.shape[0]):
        if np.sum(binary_pvalues[i,:])==1:
            cluster_assignment[i] = np.argmax(binary_pvalues[i,:])
        elif np.sum(binary_pvalues[i,:])>1:
            w = np.where(binary_pvalues[i,:])[0]
            cluster_assignment[i] = w[np.argmax(jaccard[i,w])]
        else:
            cluster_assignment[i] = np.nan

    clusters = [row_prefix+str(x) for x in np.arange(adjusted_pvalues.shape[0])]
    dictionary = dict(zip(clusters, cluster_assignment))
    dictionary = {k: v for k, v in dictionary.items() if not np.isnan(v)}
    dictionary_rows_to_cols_assignment = {v: k for k, v in dictionary.items()}
    return dictionary_rows_to_cols_assignment

def plot_densities(scores, colors_dict, names, clusters, title, figsize=(2.5, 5.5)):
    fig, axes = plt.subplots(len(clusters), 1, sharex=True, figsize=figsize)

    for ax in axes:
        pos = ax.get_position()
        # Reduce width to two-thirds (i.e. shrink by 1/3)
        ax.set_position([pos.x0, pos.y0, pos.width * 0.66, pos.height])

    for index, cluster in enumerate(clusters):

        sns.distplot(scores[cluster], color=colors_dict[cluster], kde=True,
                    hist=None, label=str(cluster), ax=axes[index], kde_kws={'linewidth': 1})
        
        # Get the line data from the KDE so you can fill between it
        l1 = axes[index].lines[0]
        x1 = l1.get_xydata()[:, 0]
        y1 = l1.get_xydata()[:, 1]
        
        axes[index].spines[['right', 'top']].set_visible(False)
        axes[index].set_title(str(cluster), loc='left', y=.3, x=.05, fontweight='bold', fontsize=8)
        axes[index].axvline(np.mean(scores[cluster]), color=colors_dict[cluster])
        axes[index].axvline(0, color='grey', linestyle='dotted')
        axes[index].fill_between(x1, y1, color=colors_dict[cluster], alpha=0.3)
        axes[index].set_ylabel('')
        
        # Process the text content
        text_str = names[cluster]
        my_list = [re.split(r' [(]GO|WP| R-', x)[0] for x in text_str]
        x_list = []
        for s in my_list:
            words = s.split()
            if len(words) > 10:
                new_str = " ".join(words[:10])
            else:
                new_str = s.strip()
            x_list.append(new_str)
        T = ('\n').join(x_list)
        
        # Add text within a colored box (placed outside the plot area)
        axes[index].text(1.15, 0.9, T,
                    transform=axes[index].transAxes,  # axes coordinates
                    fontsize=8,
                    color='black',
                    ha='left', va='top',
                    bbox=dict(facecolor='none', alpha=1, edgecolor=colors_dict[cluster]))
        
        # Add an arrow pointing from the right edge (at x=1.0) to the text box at x=1.1
        axes[index].annotate("",
                        xy=(1.1, 0.8),      # text box location
                        xycoords=axes[index].transAxes,
                        xytext=(1.0, 0.8),    # starting point at the plot edge
                        textcoords=axes[index].transAxes,
                        arrowprops=dict(arrowstyle="->", color=colors_dict[cluster], lw=1))
        
    fig.text(-0.1, 0.5, 'Density', va='center', rotation='vertical', fontsize=8)
    axes[-1].set_xlabel('gene scores [-log(p-value)*sign(logFC)]', fontsize=8)

    # Add a main title and a secondary title for the text labels
    fig.text(0.1, .9, title, va='center', fontsize=8)
    fig.text(1, .9, 'Representative cluster terms', va='center', fontsize=8)

    plt.tight_layout(rect=[0.1, 0, 1, 0.95])

