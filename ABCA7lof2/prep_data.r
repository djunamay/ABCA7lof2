sum_counts = function(counts, label, cell_labels){
    # Sums cell-level counts by factors in label vector
    #
    # Args:
    #   counts: a sparse matrix with read counts (gene x cell)
    #   label: variable of interest by which to sum counts
    #   cell_labels: vector of cell labels
    #
    # Returns:
    #   summed counts
    #   number of cells used per summation

    colnames(label) = 'ID'
    label$index = 1
    label$celltype = cell_labels
    label = as.data.frame(pivot_wider(label, values_from = index, names_from = ID))
    label[is.na(label)] = 0
    label$celltype = NULL
    summed_counts = counts%*%as.matrix(label)

    ncells = colSums(label)

    return(list('summed_counts' = summed_counts, 'ncells' = ncells))
}

get_expressed_genes = function(fraction_detected_genes_cell, fraction){
    fraction_detected_genes_cell_binary = fraction_detected_genes_cell>fraction
    genes = rownames(fraction_detected_genes_cell)
    expressed = lapply(colnames(fraction_detected_genes_cell_binary), function(x) genes[unname(fraction_detected_genes_cell_binary[,x])])
    names(expressed) = colnames(fraction_detected_genes_cell_binary)
    return(expressed)
}