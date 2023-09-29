sum_counts = function(counts, label, cell_labels){
    
  # Sums cell-level counts by factors in the label vector
  
  # Args:
  #   counts: a sparse matrix with read counts (gene x cell)
  #   label: variable of interest by which to sum counts
  #   cell_labels: vector of cell labels
  
  # Returns:
  #   A list with the following elements:
  #   - 'summed_counts': Summed counts after aggregation
  #   - 'ncells': Number of cells used per summation
    
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
    
  # Retrieves expressed genes based on a specified detection threshold
  
  # Args:
  #   fraction_detected_genes_cell: Matrix of gene detection fractions (gene x cell)
  #   fraction: Detection threshold for considering a gene as expressed
  
  # Returns:
  #   A list containing expressed genes for each cell

    fraction_detected_genes_cell_binary = fraction_detected_genes_cell>fraction
    genes = rownames(fraction_detected_genes_cell)
    expressed = lapply(colnames(fraction_detected_genes_cell_binary), function(x) genes[unname(fraction_detected_genes_cell_binary[,x])])
    names(expressed) = colnames(fraction_detected_genes_cell_binary)
    return(expressed)
}
                       
compute_stats = function(sce, sample1, sample2){
    
  # Compute statistical tests and statistics for differential gene expression analysis
  
  # Args:
  #   sce: SingleCellExperiment object containing omics expression data
  #   sample1: Sample identifier for the first group
  #   sample2: Sample identifier for the second group
  
  # Returns:
  #   A data frame with the following columns:
  #   - 'name': feature names
  #   - 'pvals': P-values from t-tests for each feature
  #   - 'logfc': Log2-fold change in expression between the two groups
  #   - 'score': Score computed as -log10(p-value) with sign based on logfc direction

    logfc = c()
    pvals = c()
    x = counts(sce)[,colData(sce)$Genotype==sample1]
    y = counts(sce)[,colData(sce)$Genotype==sample2]
    for(i in 1:nrow(x)){
        f = t.test(x[i,],y[i,])
        pvals = c(pvals, f$p.value)
        logfc = c(logfc, log2(f$estimate[[2]]/f$estimate[[1]]) )
    }
    df = as.data.frame(cbind(rownames(sce), (pvals), (logfc)))
    colnames(df) = c('name', 'pvals', 'logfc')
    df$pvals = as.numeric(df$pvals)
    df$logfc = as.numeric(df$logfc)
    df$score = sign(df$logfc)*-log10(df$pvals)
    return(df)
}

get_fatty_acid_info <- function(sce, annotation, lipid.class) {
  # Extracts fatty acid information from lipid annotation
  
  # Args:
  #   sce: SingleCellExperiment object containing expression data
  #   annotation: Name of the annotation column containing lipid information
  #   lipid.class: Name of the column indicating lipid class
  
  # Returns:
  #   A modified SingleCellExperiment object with additional metadata columns:
  #   - 'total_unsaturation': Total saturation value for each lipid
  #   - 'mean_length': Mean carbon chain length for each lipid
  #   - 'total_length': Total carbon chain length for each lipid
  #   - 'names': Names of lipids
  #   - 'merge_index': Merged index combining lipid class, total length, and total unsaturation
  #   - 'max_length': Maximum carbon chain length for each lipid
  #   - 'min_length': Minimum carbon chain length for each lipid
  #   - 'lipid.class': Lipid class for each lipid
  
  # Extract the annotation column as a data frame
  df <- as.data.frame(rowData(sce)[, c(annotation), drop = FALSE])
  
  # Initialize vectors to store calculated values
  total_saturation <- c()
  mean_length <- c()
  total_length <- c()
  max_length <- c()
  min_length <- c()
  
  # Loop through each lipid annotation
  for (i in df[[annotation]]) {
    # Split the annotation string
    temp <- strsplit(i, '_')[[1]]
    
    # Split the length and saturation components
    temp <- lapply(temp, function(x) strsplit(x, ':'))
    
    # Extract and process saturation values
    length <- unlist(lapply(temp, function(x) x[[1]][1]))
    saturation <- unlist(lapply(temp, function(x) x[[1]][2]))
    
    # Further process saturation values to remove prefixes and parentheses
    saturation <- unlist(strsplit(saturation, '[)]'))
    for (prefix in c('a', 't', 'p', 'Q', 'd', 'm', '[+O]')) {
      saturation <- (unlist(strsplit(saturation, prefix)))
    }
    
    # Further process length values to remove prefixes and parentheses
    length <- unlist(strsplit(length, '[(]'))
    for (prefix in c('a', 't', 'p', 'Q', 'd', 'm')) {
      length <- (unlist(strsplit(length, prefix)))
    }
    
    # Handle cases where length has extra parentheses
    if (length(length) == 1) {
      length <- unlist(strsplit(length, '[)]'))
    }
    
    # Convert processed values to numeric, removing NA values
    length <- (as.numeric(na.omit(length)))
    saturation <- (as.numeric(na.omit(saturation)))
    
    # Calculate and store total saturation, mean length, total length, max length, and min length
    total_saturation <- c(total_saturation, sum(saturation))
    mean_length <- c(mean_length, mean(length))
    total_length <- c(total_length, sum(length))
    max_length <- c(max_length, max(length))
    min_length <- c(min_length, min(length))
  }
  
  # Add calculated metadata columns to the original SingleCellExperiment object
  rowData(sce)$total_unsaturation <- total_saturation
  rowData(sce)$mean_length <- mean_length
  rowData(sce)$total_length <- total_length
  rowData(sce)$names <- lapply(strsplit(rownames(df), '[+]|[-]'), function(i) i[1])
  rowData(sce)$merge_index <- paste0(rowData(sce)[[lipid.class]], '_', rowData(sce)$total_length, '_', rowData(sce)$total_unsaturation)
  rowData(sce)$max_length <- max_length
  rowData(sce)$min_length <- min_length
  rowData(sce)$lipid.class <- rowData(sce)[[lipid.class]]
  
  # Return the modified SingleCellExperiment object
  return(sce)
}
                         
