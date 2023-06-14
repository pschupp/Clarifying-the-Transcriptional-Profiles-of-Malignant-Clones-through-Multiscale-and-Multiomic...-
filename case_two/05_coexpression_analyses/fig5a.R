# load in module eigengene 
# do correlation using biweight midcor as measure
# create dendrogram
# make aesthetic changes to said dendrogram
library('data.table')
library('dendextend')
# library('WGCNA')
setwd('~/@patrick/SF10711/integration_analysis/network_deconvolution/')
me = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/Module_eigengenes_02-55-42.csv')
meDist = hclust(as.dist(1-cor(me[,-1])), method='average')
meDist2 = as.dendrogram(meDist)
meDist2 = set(meDist2, 'branches_lwd', 4)
meDist2 = set(meDist2, 'leaves_pch', 21)
meDist2 = set(meDist2, 'leaves_cex', 1.5)
meDist2 = set(meDist2, 'leaves_bg', labels(meDist2))
meDist2 = set(meDist2, 'leaves_col', 'black')

pdf('dendrogram.pdf', height=3.5, width=16)
par(family='NimbusSan', font=2, cex.main=2.5, cex.lab=1.8, cex.axis=1.4, font.lab=2, lwd=3)
plot(meDist2, ylab='1 - pearson correlation', main = 'Clustering of gene coexpression modules')
axis(side = 2, lwd = 3.8) 
dev.off()
