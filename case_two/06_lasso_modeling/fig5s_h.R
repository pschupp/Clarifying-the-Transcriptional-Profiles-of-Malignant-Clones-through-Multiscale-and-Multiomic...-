library('data.table')
library('ggplot2')
dir = '/home/patrick/@patrick/SF10711/integration_analysis/lasso_modeling/glasso/cumulative/'
base::load(paste0(dir, 'workspace_after_bootstrap.Robj'))
base::load(paste0(dir, 'workspace_after_bootstrapPerm.Robj'))
plotDf=rbind(data.frame(Permuted=rep('Real', length(stability)), Stability=stability), data.frame(Permuted=rep('Permuted', length(stabilityPerm)), Stability=stabilityPerm))
plotDf$Stability=as.numeric(plotDf$Stability)
plotDf$Permuted=factor(plotDf$Permuted, levels=c('Real', 'Permuted'))
p = ggplot(plotDf, aes(x=Stability, fill=Permuted, color=Permuted)) +
    geom_density(, position="identity", alpha=.9)+
    geom_vline(aes(xintercept=45),linetype="dashed", color='red', size=2)+
    scale_fill_manual(values=c('black', 'white')) +
    scale_color_manual(values=c('black', 'black')) +
    labs(x='Stability', y='Density', title='Real and permuted stability metrics', subtitle='95th percentile cutoff of permuted is 45')+
    scale_x_continuous(limits=c(0,100), expand = c(0, 0)) +
    scale_y_continuous(limits=c(0,0.09), breaks=seq(0,0.08, 0.04), expand = c(0, 0)) +
    theme(legend.spacing.y = unit(1, 'cm')) +
    theme(legend.position="right") +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme_classic() +
    oldham_theme()

pdf('stability.plot.pdf', height=13, width=13)
p
dev.off()
