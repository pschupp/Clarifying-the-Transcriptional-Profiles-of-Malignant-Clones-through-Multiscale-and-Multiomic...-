library('data.table')
library('dendextend')
# library('WGCNA')
library('ggplot2')
library('RColorBrewer')
barOut = read.csv('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv')
setwd('~/@patrick/SF10711/integration_analysis/network_deconvolution/')
me = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/Module_eigengenes_02-55-42.csv', data.table = FALSE)
me2=fread('~/@patrick/SF10711//figures/fig6/old/alt/Bicor-None_signum0.33_minSize20_merge_ME_0.8_20246/Module_eigengenes_03-30-05.csv', data.table = FALSE)
me$lightcyan=me2$lightcyan
# trim module eigenes to sections for which we clonal abundance data (amp-seq data)
me=me[which(me$Sample %in% barOut$ampseq.sample),]
meDist = hclust(as.dist(1-cor(me[,-1])), method='average')
meDist2 = as.dendrogram(meDist)
meDist2 = set(meDist2, 'branches_lwd', 4)
meDist2 = set(meDist2, 'leaves_pch', 21)
meDist2 = set(meDist2, 'leaves_cex', 1.5)
meDist2 = set(meDist2, 'leaves_bg', labels(meDist2))
meDist2 = set(meDist2, 'leaves_col', 'black')

pdf('dendrogram.pdf', height=3.5, width=16)
par(family='NimbusSan', font=2, cex.main=2.5, cex.lab=1.8, cex.axis=1.4, font.lab=2, lwd=3)
plot(meDist2, ylab='1 - Pearson correlation', main = 'Clustering of gene coexpression modules (case 2)')
axis(side = 2, lwd = 3.8) 
dev.off()

scaleP = 2
mePlot=data.frame(Sample = me[,1],(apply(me[,-1],2,scale)))
mePlot = mePlot[, c(1, meDist$order+1)]
mePlot=reshape2::melt(mePlot, id.var='Sample')
mePlot$value[mePlot$value>2] = 2
mePlot$value[mePlot$value<(-2)] = -2
mePlot$Sample = factor(mePlot$Sample, levels = unique(mePlot$Sample))
pdf('me_heatmap.pdf', height=2.7*scaleP, width=16*scaleP)
print(ggplot(mePlot, aes(y=Sample, x=variable, fill=value))+
    geom_tile()+
#    scale_fill_distiller(palette='RdBu', limits=c(-1,1), breaks=seq(-1,1,1))+
#    scale_fill_gradient2(low='#053061', high='#67001f', mid='white')+
    scale_fill_gradientn(colors= brewer.pal(11, 'RdBu')[11:1])+
    theme_classic()+
    oldham_theme()+
    theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
 	      axis.line.x = element_blank(),
 	      axis.line.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=5),
          legend.position = 'left')+
    scale_y_discrete(expand = c(0, 0), breaks=c('SF10711_9_1_4', 'SF10711_9_2_129'), labels=c(1,140), limits=rev) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(title = 'Module eigengenes (ME)', x='', y='Section ID', fill='A.U.')
    )
dev.off()

kME=fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')

membership=data.frame(reshape2::melt(table(kME$'TopModPosFDR_0.12')))
membership$value=log10(membership$value)
membership$Var1=factor(membership$Var1, levels= meDist$labels[meDist$order])
membership=membership[match(unique(colnames(me)[-1]), membership$Var1),]
membership$yval=rep('a', nrow(membership))
pdf('me_membership.pdf', height=1.8*scaleP, width=16*scaleP)
print(ggplot(membership, aes(x=Var1, y=yval, fill=value))+
    geom_tile(color='black', size=3)+
    scale_fill_distiller(palette='OrRd', limits=c(1,4), breaks=seq(1,4), direction=1)+
    theme_classic()+
    oldham_theme()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
 	      axis.line.x = element_blank(),
 	      axis.line.y = element_blank(),
          legend.position = 'left')+
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = 'Number of genes', x='', y='', fill='Log10')
    )
dev.off()
