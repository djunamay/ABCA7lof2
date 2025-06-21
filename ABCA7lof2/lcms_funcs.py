
import os
import urllib.request
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re
from matplotlib.patches import Rectangle
def return_stats(selected_rows, names_no_blank, group1_name, group2_name, equal_var=False):
    # Create a dictionary to store t-test results for each row
    t_test_results = {}

    # Iterate over each row in the data
    for index, row in selected_rows.iterrows():
        # Get the values for the current row
        values = row.values.astype(float)
        # Separate the values into two groups based on 'names_no_blank'
        group1 = values[names_no_blank == group1_name]
        group2 = values[names_no_blank == group2_name]
        
        # Perform t-test between the two groups
        t_stat, p_value = ttest_ind(group2, group1, nan_policy='omit', equal_var=equal_var)
        
        # Calculate log2 fold change
        mean_group1 = np.mean(group1)
        mean_group2 = np.mean(group2)
        log2_fold_change = np.log2(mean_group2 / mean_group1) if mean_group1 != 0 else np.inf
        
        # Store the t-test results and log2 fold change
        t_test_results[index] = {'t_stat': t_stat, 'p_value': p_value, 'log2_fold_change': log2_fold_change}

    # Convert the results to a DataFrame for easier viewing
    t_test_results_df = pd.DataFrame.from_dict(t_test_results, orient='index')


    # Perform FDR correction on p-values
    _, fdr_corrected_pvals, _, _ = multipletests(t_test_results_df['p_value'], alpha=0.05, method='fdr_bh')

    # Add the FDR corrected p-values to the DataFrame
    t_test_results_df['fdr_corrected_p_value'] = fdr_corrected_pvals

    return t_test_results_df

def plot_volcano(pval_col, name_col, group_col,lfc_col, data, palette, force_text, max_size=50, min_size=10, alpha=0.5, subset_top_genes_by_lfc=False, label_top_genes=True, label_size=10):
    plt.figure(figsize=(5, 4))

    if subset_top_genes_by_lfc:
        top_genes = data[(data[pval_col]<0.05) & (np.abs(data[lfc_col])>1)][name_col]
    else:
        top_genes = data[data[pval_col]<0.05][name_col]

    data.loc[:,'size'] = np.where((data[pval_col]<0.05) & (np.abs(data[lfc_col])>1), max_size, min_size)
    
    sns.scatterplot(data=data, x=lfc_col, y=-1*np.log10(data[pval_col]), hue=group_col, size='size', alpha=alpha, sizes=(min_size, max_size), palette=palette, edgecolor='black')

    if label_top_genes:
        # Label top genes
        texts = []
        for gene in top_genes:
            
            if (data.loc[gene,pval_col] < 0.05): 
                color = palette[data.loc[gene, group_col]]
                txt = plt.text(data.loc[gene, lfc_col],
                            -1 * np.log10(data.loc[gene, pval_col]),
                            gene.split('+')[0], fontsize=label_size, ha='right', va='bottom', style='italic', color=color)
                texts.append(txt)

        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5), force_text=force_text)

    sns.despine(top=True, right=True)
    plt.xlabel('log$_{2}$(Fold Change)')
    plt.ylabel('-log$_{10}$(p-value) * sign(logFC)')
    plt.legend(bbox_to_anchor=(1.25, 0.5), loc='center left')
    plt.legend().remove()
    plt.title('', fontsize=10)

def plot_pca(selected_rows, names_sele, palette, draw_median=True):
    # perform PCA
    scaler = StandardScaler()
    norm_areas_subset_no_background_scaled = scaler.fit_transform(selected_rows.T)

    pca = PCA(n_components=4)
    principal_components = pca.fit_transform(norm_areas_subset_no_background_scaled)
    principal_components = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2', 'PC3', 'PC4'])
    principal_components['condition'] = names_sele

    # plot PCA
    plt.figure(figsize=(2.5, 2.5))
   
    sns.scatterplot(x='PC1', y='PC2', hue='condition', data=principal_components, palette=palette, edgecolor='w')
    explained_variance = pca.explained_variance_ratio_

    plt.xlabel(f'PC 1 ({explained_variance[0]*100:.2f}%)')
    plt.ylabel(f'PC 2 ({explained_variance[1]*100:.2f}%)')
    plt.legend(title='Line', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Remove the existing legend
    plt.legend([],[], frameon=False)
    # Add vertical line along the mean of PC1
    if draw_median:
        mean_pc1 = principal_components['PC1'].median()
        plt.axvline(x=mean_pc1, color='black', linestyle='--',  linewidth=0.5)

def plot_lfcs(heatmap_data_list, lfc_column):
    # Number of heatmaps
    n_heatmaps = len(heatmap_data_list)

    # Create a figure with n_heatmaps rows and 1 column
    fig, axes = plt.subplots(nrows=n_heatmaps, ncols=1, figsize=(2, n_heatmaps * .3))

    # If there's only one heatmap, ensure axes is iterable
    if n_heatmaps == 1:
        axes = [axes]

    # Loop over each dataframe and axis
    for ax, df in zip(axes, heatmap_data_list):
        # Sort the data by the column of interest
        sorted_data = df[[lfc_column]].sort_values(by=lfc_column)
        # Transpose so the heatmap is horizontal (one row rather than one column)
        sorted_data = sorted_data.T
        
        # Plot the horizontal heatmap on the given axis
        im = sns.heatmap(sorted_data, annot=False, cmap='RdBu_r', vmin=-1.5, vmax=1.5,
                    cbar_kws={'label': lfc_column}, ax=ax, cbar=False)
        
        # Set a title for the heatmap (displayed above each plot)
        ax.set_title('')
        ax.set_ylabel('')
        ax.set_xlabel('')
        # Optionally remove x-axis tick labels if they correspond to species names or similar unwanted labels
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(df['class'][0], rotation=0, ha='right', va='center')
        
        rect = Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                        fill=False, color="black", linewidth=1.5)
        ax.add_patch(rect)

    cbar = fig.colorbar(im.collections[0], ax=axes, orientation='horizontal', fraction=0.05, pad=0.05)
    return cbar


def count_carbons(lipid_name):
    """
    Count the total number of carbons in a lipid species name.
    The function extracts the substring within parentheses and finds all occurrences
    of a pattern where a number precedes a colon (e.g., '18:0').
    It returns the sum of these numbers.
    """
    # Extract content inside parentheses
    m = re.search(r'\((.*?)\)', lipid_name)
    if m:
        inside = m.group(1)
        # Find all numbers preceding a colon (e.g., in "d18:0", "18" is captured)
        numbers = re.findall(r'(\d+):', inside)
        if numbers:
            return sum(int(num) for num in numbers)
    # Return 0 if no match is found (or consider returning None)
    return 0


def count_unsaturations(lipid_name):
    """
    Count the total number of unsaturated bonds in a lipid species name.
    The function extracts the substring within parentheses and finds all occurrences
    of a pattern where a colon is followed by a number (e.g., '10:1').
    It returns the sum of these numbers.
    """
    # Extract the content inside parentheses
    m = re.search(r'\((.*?)\)', lipid_name)
    if m:
        inside = m.group(1)
        # Find all numbers that follow a colon (e.g., in "10:1", captures "1")
        unsat_numbers = re.findall(r':(\d+)', inside)
        if unsat_numbers:
            return sum(int(num) for num in unsat_numbers)
    # Return 0 if no unsaturation information is found
    return 0


def plot_carbon_distributions(stats):
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(5, 5))

    # TG Carbons
    axes[0, 0].hist(stats['carbons'][stats['class'] == 'TG'], bins=stats['carbons'][stats['class'] == 'TG'].nunique(), width=0.9)
    axes[0, 0].set_xlabel('Carbons')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('TG Carbons Distribution', fontsize=10)

    # TG Unsaturations
    axes[0, 1].hist(stats['unsaturations'][stats['class'] == 'TG'], bins=stats['unsaturations'][stats['class'] == 'TG'].nunique(), width=0.9)
    axes[0, 1].set_xlabel('Unsaturations')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('TG Unsaturations Distribution', fontsize=10)

    # PC Carbons
    axes[1, 0].hist(stats['carbons'][stats['class'] == 'PC'], bins=stats['carbons'][stats['class'] == 'PC'].nunique(), width=0.9)
    axes[1, 0].set_xlabel('Carbons')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('PC Carbons Distribution', fontsize=10)

    # PC Unsaturations
    axes[1, 1].hist(stats['unsaturations'][stats['class'] == 'PC'], bins=stats['unsaturations'][stats['class'] == 'PC'].nunique(), width=0.9)
    axes[1, 1].set_xlabel('Unsaturations')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title('PC Unsaturations Distribution', fontsize=10)

    # All Lipids Carbons
    axes[2, 0].hist(stats['carbons'], bins=stats['carbons'].nunique(), width=0.9)
    axes[2, 0].set_xlabel('Carbons')
    axes[2, 0].set_ylabel('Count')
    axes[2, 0].set_title('All Lipids Carbons Distribution', fontsize=10)

    # All Lipids Unsaturations
    axes[2, 1].hist(stats['unsaturations'], bins=stats['unsaturations'].nunique(), width=0.9)
    axes[2, 1].set_xlabel('Unsaturations')
    axes[2, 1].set_ylabel('Count')
    axes[2, 1].set_title('All Lipids Unsaturations Distribution', fontsize=10)

    plt.tight_layout()
    plt.show()

def plot_class_counts(stats):
    class_counts = stats['class'].value_counts()
    class_counts.plot(kind='bar')
    plt.xlabel('')
    plt.ylabel('# Detected Lipids')
    plt.title('Class Counts')
    plt.show()

def plot_with_blank(temp):
    num_plots = len(temp.columns[:-1])
    ncols = 5
    nrows = (num_plots + ncols - 1) // ncols  # Calculate the number of rows needed
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.5*ncols, 3*nrows))

    # Flatten the axes array for easy iteration
    axes = axes.flatten()

    for ax, i in zip(axes, temp.columns[:-1]):
        sns.boxplot(data=temp, x='grp', y=i, ax=ax)
        sns.swarmplot(data=temp, x='grp', y=i, ax=ax)

    # Hide any unused subplots
    for j in range(num_plots, len(axes)):
        fig.delaxes(axes[j])

    for ax in axes:
        for label in ax.get_xticklabels():
            label.set_rotation(45)
        ax.set_xlabel('')
    plt.tight_layout()
    plt.show()

def classify_fatty_acid_length(carbon_count):
    """
    Classify a fatty acid based on the number of carbons in its aliphatic tail.

    Parameters:
    carbon_count (int): Number of carbon atoms in the fatty acid chain.

    Returns:
    str: The classification of the fatty acid:
         - 'SCFA' for short-chain fatty acids (5 or fewer carbons),
         - 'MCFA' for medium-chain fatty acids (6 to 12 carbons),
         - 'LCFA' for long-chain fatty acids (13 to 21 carbons),
         - 'VLCFA' for very long-chain fatty acids (22 or more carbons).
    
    Example:
    >>> classify_fatty_acid(4)
    'SCFA'
    >>> classify_fatty_acid(8)
    'MCFA'
    >>> classify_fatty_acid(16)
    'LCFA'
    >>> classify_fatty_acid(24)
    'VLCFA'
    """
    if carbon_count <= 5:
        return 'SCFA'
    elif 6 <= carbon_count <= 12:
        return 'MCFA'
    elif 13 <= carbon_count <= 21:
        return 'LCFA'
    elif carbon_count >= 22:
        return 'VLCFA'
    else:
        return 'Unknown'

def classify_unsaturation(unsaturation):
    """
    Classify a fatty acid based on its level of unsaturation.

    Parameters:
    unsaturation (int): The number of unsaturation points (double bonds) in the fatty acid.

    Returns:
    str: The classification of the fatty acid:
         - 'UFA' if there are 0 unsaturations,
         - 'MFA' if there is 1 unsaturation,
         - 'PUFA' if there are more than 1 unsaturation.

    Example:
    >>> classify_unsaturation(0)
    'UFA'
    >>> classify_unsaturation(1)
    'MFA'
    >>> classify_unsaturation(3)
    'PUFA'
    """
    if unsaturation == 0:
        return 'SFA'
    elif unsaturation == 1:
        return 'MUFA'
    elif unsaturation > 1:
        return 'PUFA'
    else:
        return 'Unknown'