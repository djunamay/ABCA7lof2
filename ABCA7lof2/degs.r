get_gset_names_by_category = function(cat, gsets){
  # Filters gene set names by category
  
  # Args:
  #   cat: Category names to filter by
  #   gsets: List of gene set names
  
  # Returns:
  #   Filtered list of gene set names that match the specified categories
  
  gset = unlist(lapply(gsets, function(x) unlist(sum(sapply(cat, grepl, x))>0)))
  gset = (gsets[gset])
  return(gset)
}

get_limma_inputs = function(summed_counts_indexed, expressed, meta, vars){
    
  # Prepare input data for Limma analysis
  
  # Args:
  #   summed_counts_indexed: Indexed counts matrix (genes x cells) after summation
  #   expressed: List of expressed genes per cell group
  #   meta: Metadata table
  #   vars: Variables of interest from the metadata table
  
  # Returns:
  #   A list containing aggregated counts per cell group and corresponding metadata
  
    cells = strsplit(colnames(summed_counts_indexed),'[.]')
    cellnames = unlist(lapply(cells, function(x) x[[1]]))
    projids = unlist(lapply(cells, function(x) x[[2]]))
    aggs = list()
    metadata = list()
    for(i in unique(cellnames)){
      sele = projids[cellnames==i]
      curr = summed_counts_indexed[,cellnames==i]
      colnames(curr) = sele
      aggs[[i]] = curr[expressed[[i]],]
      metadata[[i]] = meta[sele,vars]
    }
    return(list('aggs' = aggs, 'metadata' = metadata))
}

RunDiffExprAnalysisLimma <- function(counts.df, var, n.sv=NULL, exclude_apoe=FALSE, exclude_batch=FALSE, exclude_both=FALSE) {
    
  # Computes differential expression using a combination of sva and voom-limma
  
  # Args:
  #   counts.df: data.frame with read counts (#genes x #samples)
  #   var: variable of interest
  #   genes.df: data.frame mapping gene IDs to gene names
  #   n.sv: number of surrogates for SVA. Automatically determined if NULL.
  #   exclude_apoe: Exclude APOE4 covariate from the model if TRUE.
  #   exclude_batch: Exclude batch effects covariate from the model if TRUE.
  #   exclude_both: Exclude both APOE4 and batch effects covariates from the model if TRUE.
  #
  # Returns:
  #   A list containing voom-limma output augmented with Storey q-values and SVA surrogate variables

    # apply edgeR normalization (TMM) to counts
    dge <- DGEList(counts=counts.df)
    dge <- calcNormFactors(dge)

    # apply voom transformation to normalized counts
    if(exclude_apoe==TRUE){
        mod1 <- model.matrix(~LOF + amyloid + nft + msex + age_death + pmi + seq_batch, data=var)
        mod0 <- model.matrix(~amyloid + nft + msex + age_death + pmi + seq_batch, data=var)
    }else if(exclude_batch==TRUE){
        mod1 <- model.matrix(~LOF + amyloid + nft + msex + age_death + pmi + APOE4, data=var)
        mod0 <- model.matrix(~amyloid + nft + msex + age_death + pmi + APOE4, data=var)
    }else if(exclude_both==TRUE){
        mod1 <- model.matrix(~LOF + amyloid + nft + msex + age_death + pmi, data=var)
        mod0 <- model.matrix(~amyloid + nft + msex + age_death + pmi, data=var)
    }else{
        mod1 <- model.matrix(~LOF + amyloid + nft + msex + age_death + pmi + seq_batch + APOE4, data=var)
        mod0 <- model.matrix(~amyloid + nft + msex + age_death + pmi + seq_batch + APOE4, data=var)
    }
    v <- voom(dge, design=mod1)

    if (is.null(n.sv)) {
        n.sv <- num.sv(v$E, mod1, method="be")
    }
    # calculate SVs on voom-transformed counts
    svobj <- sva(v$E, mod1, mod0, n.sv=n.sv)

    # re-calculate voom transformation w/ SVA covariates; run limma
    mod1 <- cbind(mod1, svobj$sv)  # intercept, var, SVs
    v <- voom(dge, design=mod1)
    fit <- lmFit(v, design=mod1)
    fit <- eBayes(fit)
    res <- topTable(fit, coef='LOF', n=Inf, sort.by="p")

    return(list("res"=res, "C"=svobj$sv))
}
      
get_deg_scores = function(degs){
  # Computes -log10(p-value) * sign(logFC) for each gene per cell type
  
  # Args:
  #   degs: output from RunDiffExprAnalysisLimma()
  #
  # Returns:
  #   A list of gene score vectors per celltype

    pseudo_scores = list()
    for(i in names(degs)){
        curr = degs[[i]]$res
        scores = sign(curr$logFC) * -log10(curr$P.Value)
        names(scores) = rownames(curr)
        scores = scores[order(scores, decreasing = TRUE)]
        pseudo_scores[[i]]$scores = scores
    }
    return(pseudo_scores)
}