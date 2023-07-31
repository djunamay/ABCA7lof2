boxplot_w_stats = function(df, x, y, group_color = x, group_fill = x, alpha=.5, palette, xlab='', ylab='', width=.5, stats_method = 'wilcox', comparisons){
    plot <- ggpubr::ggboxplot(df, x = x, y = y,
          color = group_color, fill = group_fill, alpha = alpha, palette = palette,
          xlab = xlab, ylab = ylab, width = width)+ ggpubr::stat_compare_means(method = stats_method, comparisons = comparisons) + ggplot2::geom_point(alpha = .3) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5, hjust=1))
    return(plot+scale_y_continuous(expand = expansion(mult = c(0, 0.1))))

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
