# dependencies and custom theme
# {{{
library("data.table")
library("WGCNA")
library("gplots")
library("ggplot2")
library("reshape2")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
library("scales")
library("metagMisc")
library("stringr")
library("gridExtra")
oldham_theme = function() {
    theme(
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 40, color = "black", family = "NimbusSan", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "NimbusSan", margin = margin(r = 10)),
        axis.title.y = element_text(size = 40, face = "bold", family = "NimbusSan", margin = margin(t = 0, r = -20, b = 0, l = 0)),
        axis.title.x = element_text(size = 40, face = "bold", family = "NimbusSan", margin = margin(t = -10, r = 0, b = 0, l = 0)),
        plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(t = -20, b = -40)),
        plot.subtitle = element_text(size = 40, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(t = 40, b = -40)),
        axis.line.x = element_line(size = 3),
        axis.line.y = element_line(size = 3),
        plot.margin = unit(c(4, 2, 1, 2), "lines"),
        legend.position = "right",
        legend.key.size = unit(1.3, "cm"),
        legend.text = element_text(size = 30, family = "NimbusSan")
    )
}
proc = function(muts1) {
    muts1 = muts1[, -(grep("DEL|INS", colnames(muts1)))]
    nameConv = fread("~/@patrick/SF9495/wheelplot/old.new.gene.name.conversion.txt", data.table = FALSE)
    naming = fread("~/@patrick/SF9495/wheelplot/SF9495_wheelplot_table.tsv", data.table = FALSE)
    genes = gsub("\\..*", "", colnames(muts1))[-1]
    genes = nameConv$new_gene_name[match(genes, nameConv$old_gene_name)]
    nameOut = paste(naming$Gene, gsub(".*:p\\.", "", naming$HGVSp))
    nameOut = nameOut[match(genes, naming$Gene)]
    nameOut[19] = "PHF8 3' UTR"
    nameOut = gsub("Leu", "L", (nameOut))
    nameOut = gsub("Ile", "I", (nameOut))
    nameOut = gsub("Pro", "P", (nameOut))
    nameOut = gsub("Phe", "F", (nameOut))
    nameOut = gsub("Arg", "R", (nameOut))
    nameOut = gsub("His", "H", (nameOut))
    nameOut = gsub("Ser", "S", (nameOut))
    nameOut = gsub("Cys", "C", (nameOut))
    nameOut = gsub("Val", "V", (nameOut))
    nameOut = gsub("Thr", "T", (nameOut))
    nameOut = gsub("Ala", "A", (nameOut))
    nameOut = gsub("Asp", "D", (nameOut))
    nameOut = gsub("Asn", "N", (nameOut))
    nameOut = gsub("Gly", "G", (nameOut))
    nameOut = gsub("Gln", "N", (nameOut))
    nameOut = gsub("Glu", "E", (nameOut))
    colnames(muts1)[2:ncol(muts1)] = nameOut
    muts1 = data.frame(muts1)[-nrow(muts1), ] # remove blood
    return(muts1)
}
setwd("~/@patrick/SF9495/figures")
# }}}
# read in data, tidy up text, cluster, and prepare color vectors
# {{{
muts1 = fread("~/@patrick/SF9495/ampseq/Combined_variant_frequencies_run1_and_run2_with_readCounts.csv", data.table = FALSE) 
muts1_depth = muts1[, c(1, grep("depth$", colnames(muts1)))]
muts1_freq = muts1[, c(1, grep("_freq$", colnames(muts1)))]
muts1_SE = muts1[, c(1, grep("_SE$", colnames(muts1)))]
freqs.trim = proc(muts1_freq)
depth = proc(muts1_depth)
SE = proc(muts1_SE)
freqs.trim = freqs.trim[, -grep("PHF8|PCLO", colnames(freqs.trim))]
depth = depth[, -grep("PHF8|PCLO", colnames(depth))]
SE = SE[, -grep("PHF8|PCLO", colnames(SE))]
freqs = freqs.trim
colnames(freqs)[-1] = paste0(gsub("\\.", " ", colnames(freqs)[-1]), " VAF")
colnames(depth)[-1] = paste0(gsub("\\.", " ", colnames(depth)[-1]), " depth")
colnames(SE)[-1] = paste0(gsub("\\.", " ", colnames(SE)[-1]), " SE")
table = data.table(freqs[, -1], SE[, -1], depth[, -1])
tabOrder = names(sort(apply(freqs[, -1], 2, sum), decreasing = T))
tabOrder = match(gsub("\\s.*", "", tabOrder), gsub("\\s.*", "", colnames(table)))
table.t = data.table(Section = freqs[, 1], table[, as.numeric(unlist(sapply(tabOrder, function(x) (c(x, x + 22, x + 44)))))])
write.table(table.t, file = "ampseq_table.csv", row.names = F, quote = F, sep = ",")
cluster2 = hclust(as.dist(1 - cor(freqs.trim[, -1])))
clustDist = as.dist(1 - cor(freqs.trim[, -1]))
corMat2 = cor(freqs.trim[, -1])
for (i in seq(1, nrow(corMat2))) {
    corMat2[i, i] = 1
}
groups1 = cutree(cluster2, k = 4)
groups1 = as.factor(groups1)
names(groups1) = gsub("\\.", " ", names(groups1))
cluster.numb = length(unique(groups1))
Clusters = data.frame(paste("c", groups1, sep = ""))
rownames(Clusters) = names(groups1)
colnames(Clusters) = "Clusters"
map = list(Clusters = setNames(brewer.pal(cluster.numb, "Set1"), as.character(unique(Clusters$Clusters))))
colnames(corMat2) = rownames(corMat2) = names(groups1)
corMat2 = corMat2[match(rownames(Clusters), rownames(corMat2)), match(rownames(Clusters), colnames(corMat2))]
colnames(depth) = gsub("\\.", " ", colnames(depth))
depthVafs = apply(depth[, names(groups1)[which(groups1 == 1)]], 2, mean)
mean(depthVafs[1:8])
mean(depthVafs[9:15])
breaksList = seq(-1, 1, by = 0.01)
colnames(corMat2) = gsub("\\.", " ", colnames(corMat2))
# }}}
# create heatmap
# {{{
colorText = rep("black", 22)
faceText = rep("plain", 22)
faceText[grep("ACCS|GDPD1|TP53|CTNNA3|IDH1|TEC|FAM193B|IAPP|GPR173", colnames(corMat2))] = "bold"
library("dendextend")
colorVec = brewer_pal(palette = "Set1")(9)[c(4, 3, 5)]
colFun = colorRamp2(seq(-1, 1, length = 7), rev(brewer.pal(n = 7, name = "RdBu")))
row_dend = color_branches(cluster2, k = 3, col = colorVec)
pdf("high.vaf.heatmap.dist.pearson.pdf", width = 8, height = 8)
ht = Heatmap(corMat2,
    col = colFun, rect_gp = gpar(lwd = 3, col = "black"), show_row_dend = F, column_labels = rep("", ncol(corMat2)), column_dend_height = unit(.75, "cm"),
    cluster_columns = row_dend, column_split = 3, row_order = cluster2$order, column_gap = unit(0, "mm"), column_dend_gp = gpar(lwd = 5),
    clustering_distance_columns = clustDist, show_row_names = F, column_title = NULL,
    heatmap_legend_param = list(col_fun = colFun, title = "Pearson correlation", direction = "horizontal", at = seq(-1, 1, 1), border = F, legend_width = unit(8, "cm"), title_gp = gpar(fontsize = 20, fontface = "bold", fontfamily = "NimbusSan"), title_position = "topcenter", labels_gp = gpar(fontsize = 20, fontfamily = "NimbusSan")),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if ((i == 12) & (j == 12)) grid.text("3", x, y, gp = gpar(fontsize = 40, fontfamily = "NimbusSan", fontface = "bold", col = "white"))
        if ((i == 2) & (j == 2)) grid.text("1", x, y, gp = gpar(fontsize = 40, fontfamily = "NimbusSan", fontface = "bold", col = "white"))
        if ((i == 11) & (j == 11)) grid.text("2", x, y, gp = gpar(fontsize = 40, fontfamily = "NimbusSan", fontface = "bold", col = "white"))
    },
    top_annotation = HeatmapAnnotation(ta = anno_block(gp = gpar(fill = colorVec, lwd = 5), labels = c("Cl. 2", "Cl. 3", "Cl. 1"), labels_gp = gpar(fontsize = 20, fontface = "bold", fontfamily = "NimbusSan")))
) +
    rowAnnotation(ra = anno_text(rownames(corMat2), just = "left", gp = gpar(col = colorText, fontface = faceText, fontfamily = "NimbusSan", fontsize = 14)))
draw(ht, heatmap_legend_side = "bottom", column_title = "Clustering of Amp-seq VAFs", column_title_gp = gpar(fontsize = 28, fontface = "bold", fontfamily = "NimbusSan"))
dev.off()
# }}}
# prepare data for lineplot by creating the appropriate clusters and adding mean depth
# {{{
colnames(freqs.trim) = gsub("\\.", " ", colnames(freqs.trim))
groups1.bak = groups1
groups1 = as.numeric(groups1)
names(groups1) = names(groups1.bak)
groups1[grep("ACCS|GDPD1|TP53|CTNNA3|IDH1|TEC|FAM193B|IAPP", names(groups1))] = 0
groups1[grep("SUSD4", names(groups1))] = 1
titles = c("Cluster 1", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
sub = c(expression(paste("mean depth: 1.6x10"^"4 ", "reads")), expression(paste("mean depth: 3.7x10"^"3 ", "reads")), expression(paste("mean depth: 4.4x10"^"3 ", "reads")), expression(paste("mean depth: 2.9x10"^"3 ", "reads")), "")
freqs.trim = freqs.trim[, c(1, cluster2$order + 1)]
# }}}
# create lineplot
# {{{
lims = list(c(0, .9), c(0, .9), c(0, .3), c(0, .2))
pdf("lineplot_clusters.pdf", width = 24, height = 19)
i = 1
p1 = list()
p1[[1]] = 1
p1[[2]] = 2
p1[[3]] = 3
p1[[4]] = 4
for (ea in seq(0, 3)) {
    freqs.clust1 = freqs.trim[, c(1, which(colnames(freqs.trim) %in% names(groups1[which(groups1 == ea)])))]
    if (ncol(freqs.clust1) == 2) {
        freqs.clust1.m = data.frame(sample = freqs.clust1$sample, variable = rep(colnames(freqs.clust1)[2], nrow(freqs.clust1)), value = freqs.clust1[, 2])
    } else {
        freqs.clust1.m = melt(freqs.clust1, id.vars = "sample")
    }
    freqs.clust1.m[, 1] = as.numeric(freqs.clust1.m[, 1])
    freqs.clust1.m[, 2] = as.factor(freqs.clust1.m[, 2])
    freqs.clust1.m[, 2] = factor(freqs.clust1.m[, 2], levels = as.character(unique(freqs.clust1.m[, 2])))
    if (i == 1) {
        p1[[i]] = (ggplot(freqs.clust1.m, aes(x = sample, y = value, group = variable, color = variable)) +
            theme_classic() +
            oldham_theme() +
            theme(plot.title = element_text(hjust = 1), plot.subtitle = element_text(hjust = -2.2)) +
            geom_line(size = 2) +
            labs(x = "Section ID", y = "VAF", title = titles[i], subtitle = sub[i]) +
            scale_x_continuous(expand = c(0, 0), limits = c(1, 81), breaks = c(1, 81), labels = c(1, 81)) +
            scale_y_continuous(expand = c(0, 0), limits = lims[[1]], breaks = lims[[1]], labels = lims[[1]]) +
            theme(legend.text = element_text(face = "bold"), plot.title = element_text(color = brewer_pal(palette = "Set1")(9)[5])) +
            scale_color_brewer(palette = "Set1"))
    } else {
        p1[[i]] = (ggplot(freqs.clust1.m, aes(x = sample, y = value, group = variable, color = variable)) +
            theme_classic() +
            oldham_theme() +
            theme(plot.title = element_text(hjust = 1), plot.subtitle = element_text(hjust = -2.2)) +
            geom_line(size = 2) +
            labs(x = "Section ID", y = "VAF", title = titles[i], subtitle = sub[i]) +
            scale_x_continuous(expand = c(0, 0), limits = c(1, 81), breaks = c(1, 81), labels = c(1, 81)) +
            scale_y_continuous(expand = c(0, 0), limits = lims[[i]], breaks = lims[[i]], labels = lims[[i]]) +
            theme(plot.title = element_text(color = colorVec[c(3, 1, 2)][i - 1])) +
            scale_color_brewer(palette = "Set1"))
    }
    i = i + 1
}
grid.arrange(p1[[1]] + theme(plot.margin = unit(c(1, 1, 0, 0.5), "cm")), p1[[2]] + theme(plot.margin = unit(c(1, 1, 0, 0.5), "cm")), p1[[3]] + theme(plot.margin = unit(c(1, 1, 0, 0.5), "cm")), p1[[4]] + theme(plot.margin = unit(c(1, 1, 0, 0.5), "cm")), ncol = 2, padding = unit(0, "line"))
dev.off()
# }}}
