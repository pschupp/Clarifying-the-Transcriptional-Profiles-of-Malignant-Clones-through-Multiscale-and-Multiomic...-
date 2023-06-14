# load dependencies
# {{{
library("ggplot2")
library("data.table")
library("dendextend")
library("WGCNA")
library("data.table")
library("ggplot2")
library("RColorBrewer")
WD = "~/@patrick/SF9495/figures/"
# }}}
# create dendrogram
# {{{
me = fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/Module_eigengenes_06-26-41.csv", data.table = FALSE)
meDist = hclust(as.dist(1 - cor(me[, -1])), method = "average")
meDist2 = as.dendrogram(meDist)
meDist2 = set(meDist2, "branches_lwd", 4)
meDist2 = set(meDist2, "leaves_pch", 21)
meDist2 = set(meDist2, "leaves_cex", 1.5)
meDist2 = set(meDist2, "leaves_bg", labels(meDist2))
meDist2 = set(meDist2, "leaves_col", "black")

pdf(paste0(WD, "dendrogram.pdf"), height = 3.5, width = 16)
par(family = "NimbusSan", font = 2, cex.main = 2.5, cex.lab = 1.8, cex.axis = 1.4, font.lab = 2, lwd = 3)
plot(meDist2, ylab = "1 - Pearson correlation", main = "Clustering of gene coexpression modules (case 1)")
axis(side = 2, lwd = 3.8)
dev.off()
# }}}
# create heatmap of module eigengenes
# {{{
scaleP = 2
me = fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/Module_eigengenes_06-26-41.csv", data.table = FALSE)[-c(68, 69), ]
mePlot = data.frame(Sample = me[, 1], (apply(me[, -1], 2, scale)))
mePlot = reshape2::melt(mePlot, id.var = "Sample")
mePlot$value[mePlot$value > 2] = 2
mePlot$value[mePlot$value < (-2)] = -2
pdf(paste0(WD, "me_heatmap.pdf"), height = 2.7 * scaleP, width = 16 * scaleP)
print(ggplot(mePlot, aes(y = Sample, x = variable, fill = value)) +
    geom_tile() +
    #    scale_fill_distiller(palette='RdBu', limits=c(-1,1), breaks=seq(-1,1,1))+
    #    scale_fill_gradient2(low='#053061', high='#67001f', mid='white')+
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu")[11:1]) +
    theme_classic() +
    oldham_theme() +
    theme(
        axis.text.x = element_blank(),
        #         axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 5),
        legend.position = "left"
    ) +
    scale_y_discrete(expand = c(0, 0), breaks = c("SF9495_1", "SF9495_67"), labels = c(1, 81), limits = rev) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(title = "Module eigengenes (ME)", x = "", y = "Section ID", fill = "A.U."))
dev.off()
# }}} 
# create heatmap of number of members per module
# {{{
kME = fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/kME_table_06-26-41.csv")
scaleP = 2
# membership=data.frame(reshape2::melt(table(kME$ModSeed)))
membership = data.frame(reshape2::melt(table(kME$TopModPosFDR_0.0549)))
membership$value = log10(membership$value)
membership$Var1 = factor(membership$Var1, levels = unique(colnames(me)[-1]))
membership = membership[match(unique(colnames(me)[-1]), membership$Var1), ]
membership$yval = rep("a", nrow(membership))
pdf(paste0(WD, "me_membership.pdf"), height = 1.8 * scaleP, width = 16 * scaleP)
print(ggplot(membership, aes(x = Var1, y = yval, fill = value)) +
    geom_tile(color = "black", size = 3) +
    scale_fill_distiller(palette = "OrRd", limits = c(1, 4), breaks = seq(1, 4), direction = 1) +
    theme_classic() +
    oldham_theme() +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "left"
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = "Number of genes", x = "", y = "", fill = "Log10"))
dev.off()
# }}}
