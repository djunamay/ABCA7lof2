boxplot_w_stats = function(df, x, y, group_color = x, group_fill = x, alpha=.5, palette, xlab='', ylab='', width=.5, stats_method = 'wilcox', comparisons){
    plot <- ggpubr::ggboxplot(df, x = x, y = y,
          color = group_color, fill = group_fill, alpha = alpha, palette = palette,
          xlab = xlab, ylab = ylab, width = width)+ ggpubr::stat_compare_means(method = stats_method, comparisons = comparisons) + ggplot2::geom_point(alpha = .3) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5, hjust=1))
    return(plot+scale_y_continuous(expand = expansion(mult = c(0, 0.1))))

}

plot_coords_by_gene = function (df, x_name, y_name, annotation_name, alpha_name, low, high, s=.0001) 
{
    df = (df[,c(x_name, y_name, annotation_name, alpha_name)])
    colnames(df) = c('x_name', 'y_name', 'annotation_name','alpha_name')
    df$x_name = as.numeric(df$x_name)
    df$y_name = as.numeric(df$y_name)
    df$annotation_name = as.numeric(df$annotation_name)
    df$alpha_name = as.numeric(df$alpha_name)
    
    plot = ggplot(df, aes(x=x_name, y=y_name)) + 
                  geom_point(size = s, aes( alpha = alpha_name, color = annotation_name), show.legend=TRUE) + theme_void()+theme( panel.background = element_rect(colour = "black", size=0), legend.position="none") +  
                  scale_color_gradient(low=low, high = high ) + theme(text = element_text(size=10), legend.position="right")  + guides(alpha=FALSE) + labs(color="expression")
    return(plot)
}

get_permutation_plot = function(gene, df){
    rbfox_null = c()

    for(i in 1:10000){   
        x = sample(df$LOF)
        rbfox_null = c(rbfox_null, mean(df[,gene][x=='LoF'])-mean(df[,gene][x=='Con']))
    }
    rbfox_null = rbfox_null[order(abs(rbfox_null), decreasing = T)]
    get_pval = approxfun(x = rbfox_null, y = (1:length(rbfox_null))/length(rbfox_null))

    observation = mean(df[,gene][df$LOF=='LoF'])-mean(df[,gene][df$LOF=='Con'])
    p <- ggplot(as.data.frame(rbfox_null), aes(x=rbfox_null)) + 
      geom_histogram(bins = 100) + theme_classic() + geom_vline(xintercept =observation, color="red", linetype="dashed", size=1) + annotate(geom="text", x=.1, y=300, label=paste0(' ',round(get_pval(observation), digits = 3)),
                  color="red") + ggtitle(gene) + ylab('counts') + xlab('E(ABCA7 LoF) - E(ABCA7 Con)')
    return(p)
}


get_barplot = function(all_data, x, y){

    df = as.data.frame(all_data[,c(x,y)])
    colnames(df) = c('var1', 'var2')

    res = stats::fisher.test(df$var1, df$var2)
    txt = paste0('  ',round(res$p.value, digits = 2))

    bp = ggplot(df) +
      aes(x = factor(var1), fill = factor(var2)) +
      geom_bar(color = 'black', position = "fill") + ggtitle(y) + annotate("text", x=levels(factor(df$var1))[1], y=1.1, label= txt) + xlab('LoF variants') + ylab('fraction') + theme_classic() + labs(x = "", fill = y) +scale_fill_brewer(palette="Dark2")

    return(bp)

}

plot_coords_by_grp = function (df, x_name, y_name, annotation_name, alpha_name, colors) 
{
    df = (df[,c(x_name, y_name, annotation_name, alpha_name)])
    colnames(df) = c('x_name', 'y_name', 'annotation_name', 'alpha_name')
    
    plot = ggplot(df, aes(x = x_name, y = y_name, color = annotation_name, 
        alpha = alpha_name)) + geom_point(size = 0.001) + theme_void() + 
        theme(panel.background = element_rect(colour = "black", 
            size = 0), legend.position = "none") + scale_color_manual(values = (colors[as.character(unique(df$annotation_name))])) + 
        ggtitle("") + theme(text = element_text(size = (11))) + 
        theme(legend.position = "right") + guides(colour = guide_legend(override.aes = list(size = 3))) + 
        theme(legend.title = element_blank(), legend.position = "none")
    return(plot + theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", 
            color = NA), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(fill = "transparent")))
    return(plot)
}


plot_degs = function(data, gene, name){
    
    data$score.x = sign(data$logFC.x) * -log10(data$P.Value.x)
    data$score.y = sign(data$logFC.y) * -log10(data$P.Value.y)
    
    x = cor.test(data$score.x, data$score.y)

    estimate = x$estimate
    lower = x$conf.int[1]
    upper = x$conf.int[2]

    text = paste0('pcc = ',round(estimate, 2), ', 95% CI: [', round(lower, 2), ', ', round(upper, 2), ']')

    ggplot(data, aes(x=as.numeric(score.x), y=as.numeric(score.y)), label = Row.names) +
    geom_hex(bins = 100) + geom_point(data = data[data$Row.names%in%gene,], col = 'red') + geom_smooth(method=lm, formula='y ~ x') + theme_classic() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)  + scale_fill_viridis_c() +
    geom_text_repel(label = ifelse(data$Row.names%in%gene,data$Row.names,''), max.overlaps  = 1000000) + xlab('Score: All Samples') + ylab('Score: Subset Samples') + ggtitle(paste0('comparison of DEGs for ', name)) +
    annotate("text", x=0, y=max(data$score.y), label= text)
}

plot_lipid_diff = function(data, gene, name, x_label, y_label){

    x = cor.test(data$sample1, data$sample2)

    estimate = x$estimate
    lower = x$conf.int[1]
    upper = x$conf.int[2]

    text = paste0('pcc = ',round(estimate, 2), ',\n 95% CI: [', round(lower, 2), ', ', round(upper, 2), ']')

    ggplot(data, aes(x=as.numeric(sample1), y=as.numeric(sample2)), label = Row.names) +
    geom_point(data = data[!(data$Row.names%in%gene),], col = 'grey') + geom_point(data = data[data$Row.names%in%gene,], col = 'forestgreen') + geom_smooth(method=lm, formula='y ~ x') + theme_classic() + geom_vline(xintercept = 0, linetype='dashed') + geom_hline(yintercept = 0, linetype='dashed')  + scale_fill_viridis_c() +
    geom_text_repel(label = ifelse(data$Row.names%in%gene,data$Row.names,''), max.overlaps  = 1000000) + xlab(x_label) + ylab(y_label) + ggtitle(name) +
    annotate("text", x=-1.5, y=max(data$sample2), label= text)
}

draw_deg_hmap = function(celltype, all_data){
    d = all_data$av_logcounts_by_ind[[celltype]][genes,] %>% melt(.) %>% as_tibble(.)
    genes = rownames(degs$degs_all[[celltype]]$res[order(degs$degs_all[[celltype]]$res$P.Value,decreasing=F),])[1:30]
    d$batch = all_data$summary[as.character(d$Var2),'seq_batch']
    d$lof = all_data$summary[as.character(d$Var2),'LOF']

    h1 = d %>%
    tidyHeatmap::heatmap(.value=value, .row=Var1, .column=Var2, scale='row', clustering_distance_columns='pearson',clustering_distance_rows='pearson', palette_value = circlize::colorRamp2(seq(-2, 2, length.out = 11), rev(RColorBrewer::brewer.pal(11, "RdBu"))), column_title=celltype,row_title='') %>%
    add_tile(lof)%>%wrap_heatmap()
    return(h1)
}

draw_lipid_heatmap = function(d, N, lfc_cut, pval_cut, p, mat, show_heatmap_legend, i){
    
    d1 = d%>%group_by(class, variable)%>%top_n(-N, score)%>%filter(logfc< -1*lfc_cut, pvals<pval_cut)%>%ungroup#%>%filter(class%in%c('Neutral lipids', 'Phospholipids', 'Sphingolipids'))
    d2 = d%>%group_by(class, variable)%>%top_n(N, score)%>%filter(logfc> lfc_cut, pvals<pval_cut)%>%ungroup#%>%filter(class%in%c('Neutral lipids', 'Phospholipids', 'Sphingolipids'))
    d = rbind(d1,d2)%>%arrange(class)
    #d$class = factor(d$class, levels = names(class_cols))
    print(length(unique(d$name)))
    cluster_cols = hclust(dist(scale(t(mat[unique(d$name),]))))
    #cluster_rows = hclust(dist(t(scale(t(mat[unique(d$name),]+0.0001)))))
    
    if (unname(unlist(unique(d[d$variable==cluster_cols$labels[1],'genotype'])))!=i){
        cluster_cols = rev(cluster_cols)
    }

    #if (unname(unlist(unique(d[d$name==cluster_rows$labels[1],'logfc'])))>0){
     #   cluster_rows = rev(cluster_rows)
    #}
    
    h1 = d %>% group_by(class)%>%
    tidyHeatmap::heatmap(
        column_title='',
        row_title='',
     .row = lipid_name,
     .column = variable,
     .value = value,
     scale = "row",show_heatmap_legend =show_heatmap_legend, show_row_dend = FALSE, show_annotation_name=FALSE,
    palette_value = circlize::colorRamp2(
            seq(-3, 3, length.out = 11), 
            rev(RColorBrewer::brewer.pal(11, "RdBu"))), 
    #width = length(unique(d$variable))*unit(5, "mm"), height = length(unique(d$lipid_name))*unit(5, "mm"), 
        cluster_columns=cluster_cols) %>%
    add_tile(class, palette = p, show_legend = show_heatmap_legend)%>%
    add_tile(genotype, palette = c('grey', 'red'), show_legend = FALSE)%>%wrap_heatmap()

    return(h1)
}


plot_volcano = function(sce, ratio_name, pval_name, lab, pval_cut, lfc_cut, adjust, annotation1, class_cols){
    cols = c('blue', 'grey', 'red')
    names(cols) = c('down', 'other', 'up')
    temp = as.data.frame(rowData(sce)[,c(ratio_name, pval_name, annotation1)])
    colnames(temp) = c('log2', 'pvalue', 'label')
    if(adjust){
        temp$pvalue = p.adjust(temp$pvalue, 'fdr')
    }
    temp$direction = ifelse(temp$log2>lfc_cut & temp$pvalue<pval_cut, 'up', ifelse(temp$log2< -1*lfc_cut & temp$pvalue<pval_cut, 'down', 'other'))
    temp$label2 = ifelse(temp$direction!='other' & temp$label%in%lab, temp$label,'')
    p = ggplot(temp, aes(x=log2, y= -log10(pvalue), col=label, label = label2))+ geom_text_repel(max.overlaps = Inf, size = 5)+  scale_color_manual(values = c('yellow', 'blue','green')) + geom_point(aes(alpha = 0.01), size = 2, show.legend = FALSE)+geom_hline(yintercept = -log10(pval_cut)) +geom_vline(xintercept = 0)+ theme_classic() +   theme(text = element_text(size=20),strip.background = element_blank(),strip.placement = "outside" ,strip.text.y = element_text(size = 20, color = "black"))+ theme(legend.position='none')#+ #facet_wrap(~label, ncol = 7, scales = "free_x")
    return(p)
}