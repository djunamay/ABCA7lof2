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
                       
compute_stats = function(sce, sample1, sample2){
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

# extract total carbon chain length and unsaturation values
get_fatty_acid_info = function(sce, annotation, lipid.class){
    df = as.data.frame(rowData(sce)[,c(annotation), drop=FALSE])
    total_saturation = c()
    mean_length = c()
    total_length = c()
    max_length = c()
    min_length = c()
    for(i in df[[annotation]]){
        temp = strsplit(i, '_')[[1]]
        temp = lapply(temp, function(x) strsplit(x, ':'))
        length = unlist(lapply(temp, function(x) x[[1]][1]))
        saturation = unlist(lapply(temp, function(x) x[[1]][2]))
        saturation = ((unlist(strsplit(saturation, '[)]'))))
        for(prefix in c('a','t','p','Q', 'd', 'm', '[+O]')){
            saturation=(unlist(strsplit(saturation, prefix)))
        }
        length = ((unlist(strsplit(length, '[(]'))))                          
        for(prefix in c('a','t','p','Q', 'd', 'm')){
            length=(unlist(strsplit(length, prefix)))
        }
        if(length(length)==1){
            length=unlist(strsplit(length, '[)]'))
        }
        length = (as.numeric(na.omit(length)))
        saturation = (as.numeric(na.omit(saturation)))
        #print(length)                        
        total_saturation = c(total_saturation, sum(saturation))
        mean_length = c(mean_length, mean(length))
        total_length = c(total_length, sum(length))
        max_length = c(max_length, max(length)) 
        min_length = c(min_length,min(length))                            
    }
    rowData(sce)$total_unsaturation = total_saturation
    rowData(sce)$mean_length = mean_length
    rowData(sce)$total_length = total_length
    rowData(sce)$names = lapply(strsplit(rownames(df), '[+]|[-]'), function(i) i[1])
    rowData(sce)$merge_index = paste0(rowData(sce)[[lipid.class]], '_', rowData(sce)$total_length, '_', rowData(sce)$total_unsaturation)#, '_', rowData(sce)$total_length)
    #df$direction = ifelse(df[[ratio_name]]>lfc_cut & df[[pval_name]]<pval_cut, 'up', ifelse(df[[ratio_name]]< -1*lfc_cut & df[[pval_name]]<pval_cut, 'down', 'other'))                               
    rowData(sce)$max_length = max_length
    rowData(sce)$min_length = min_length 
    rowData(sce)$lipid.class = rowData(sce)[[lipid.class]]
    return(sce)
}
                                
