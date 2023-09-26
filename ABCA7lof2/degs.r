get_gset_names_by_category = function(cat, gsets){
  gset = unlist(lapply(gsets, function(x) unlist(sum(sapply(cat, grepl, x))>0)))
  gset = (gsets[gset])
  return(gset)
}
                       
#################################
get_degs = function(nebula_out, p_cut, LFC_cut, p_slot, LFC_slot){
    degs = list()
    for(i in names(nebula_out)){
        df = (nebula_out[[i]]$all)
        df_up = df[df[[p_slot]]<p_cut & df[[LFC_slot]]>LFC_cut,'gene']
        df_down = df[df[[p_slot]]<p_cut & df[[LFC_slot]]<LFC_cut,'gene']
        degs[[i]][['up']] = df_up
        degs[[i]][['down']] = df_down
    }
    return(degs)
}
#################################
filter_genes_nebula = function(nebula_in_stem_path, celltype_vector, conv=-20, dispersion=1, subject_only = TRUE){

    all = list()
    for(i in celltype_vector){
        print(i)
        path = paste0(nebula_in_stem_path, '_', i, '.rds')
        temp = readRDS(path)
        keep = temp$convergence>conv
        keep2 = temp$overdispersion<dispersion
        print('after conv filter:')
        print(table(keep))
        if(subject_only == TRUE){
            print('after overdispersion filter')
            print(table(keep2))
        }else{
            keep3 = temp$overdispersion[,'Cell']<100
            keep2 = keep2 + keep3

        }
        keep = keep+keep2
        print('after both filters:')
        print(table(keep))
        all[[i]] =  temp$summary$gene[keep==2]
    }
    return(all)
}
#################################
get_scores_from_nebula = function(nebula_out){
    scores = lapply(names(nebula_out), function(x) cbind(nebula_out[[x]]$all[,c('score', 'gene')], x))
    df = do.call('rbind', scores)
    scores = as.data.frame(tidyr::pivot_wider(df, values_from = 'score', names_from = 'x', id_cols = 'gene'))
    rownames(scores) = scores$gene
    scores$gene = NULL
    scores[is.na(scores)] = 0
    return(scores)
}
#################################
get_degs = function(nebula_out, p_cut, LFC_cut, p_slot, LFC_slot){
    degs = list()
    for(i in names(nebula_out)){
        df = (nebula_out[[i]]$all)
        df_up = df[df[[p_slot]]<p_cut & df[[LFC_slot]]>LFC_cut,'gene']
        df_down = df[df[[p_slot]]<p_cut & df[[LFC_slot]]<LFC_cut,'gene']
        degs[[i]][['up']] = df_up
        degs[[i]][['down']] = df_down
    }
    return(degs)
}
#################################
RunDiffExprAnalysisLimma_pseudobulk <- function(logcounts, var, predict) {
    mod1 = model.matrix(~LOF + amyloid + tangles + msex + age_death + pmi + seq_batch + APOE4, data=predict)
    fit <- lmFit(logcounts, design=mod1)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=var, n=Inf, sort.by="none")
    res = res[!duplicated(res$ID), ]
    rownames(res) = res$ID
    return(list("res"=res))
}
#################################
get_scores = function(nebula_in_stem_path, celltype_vector, p.slot, lfc.slot, direction = TRUE){

    all = list()
    for(i in celltype_vector){

        path = paste0(nebula_in_stem_path, '_', i, '.rds')
        out = readRDS(path)
        out = out$summary
        out$padj = stats::p.adjust(out[,p.slot], method = 'fdr')
        out = out[, c('gene', p.slot, lfc.slot, 'padj')]
        out = out[!duplicated(out$gene),]
        rownames(out) = out$gene
        if(isTRUE(direction)){
            out$score = -log10(out[,p.slot])*sign(out[,lfc.slot])
            out = out[order(out$score,decreasing=TRUE),]
            all[[i]]$all = out
            x = out$score
            names(x) = out$gene
            all[[i]]$scores = x
        }else{
            out$score = -log10(out[,p.slot])
            out = out[order(out$score,decreasing=TRUE),]
            all[[i]]$all = out
            x = out$score
            names(x) = out$gene
            all[[i]]$scores = x
        }

    }

    return(all)
}
#################################


RunDiffExprAnalysisLimma <- function(counts.df, var, n.sv=NULL, exclude_apoe=FALSE, exclude_batch=FALSE, exclude_both=FALSE) {
    # Computes differential expression using a combination of sva and voom-limma
    #
    # Args:
    #   counts.df: data.frame with read counts (#genes x #samples)
    #   var: variable of interest
    #   genes.df: data.frame mapping gene IDs to gene names
    #   n.sv: number of surrogates for SVA. Automatically determined if NULL.
    #
    # Returns:
    #   voom-limma output augmented with Storey q-values
    #   SVA surrogate variables


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

get_limma_inputs = function(summed_counts_indexed, expressed, meta, vars){
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

#################################
get_deg_scores = function(degs){
    # computes -log10(p-value) * sign(logFC) for each gene per cell type
    #
    # Args:
    #   degs: output from RunDiffExprAnalysisLimma()
    #
    # Returns:
    #   list of gene score vectors per celltype

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