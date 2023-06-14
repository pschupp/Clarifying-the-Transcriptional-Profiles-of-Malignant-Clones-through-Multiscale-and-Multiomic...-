library("data.table")
library("ggplot2")
# read in cellular abundance
barOut = fread("~/@patrick/SF9495/figures/sf9495_cellular_abundance.csv")
colnames(barOutIn) = c("clone.1", "clone.2", "clone.3", "clone.4", "clone.5", "clone.6")
barOutIn = data.frame(barOut[, -c(1, 2, 3, 4, 5)], Malignant = 1 - barOut$Nonmalignant, Clone.2 = (barOut$clone.2 + barOut$clone.3 + barOut$clone.4 + barOut$clone.5 + barOut$clone.6), Clone.3 = (barOut$clone.3 + barOut$clone.4 + barOut$clone.5))
colnames(barOutIn)[1:3] = c("Clone.4", "Clone.5", "Clone.6")
barOutIn = barOutIn[, c(4, 5, 6, 1, 2, 3)]
colnames(barOutIn) = c("Clone.1", "Clone.2", "Clone.3", "Clone.4", "Clone.5", "Clone.6")
barOutIn = barOutIn[, -2]
corPlot = reshape2::melt(cor(barOutIn))
corPlot = data.table(corPlot, label = round(corPlot$value, 2))
corPlot$Var1 = gsub("\\.", " ", corPlot$Var1)
corPlot$Var2 = gsub("\\.", " ", corPlot$Var2)
setwd("~/@patrick/SF9495/integration_analysis/network_deconvolution")
# create correlation heatmap, cumulative version
pdf("cor_heatmap.pdf", height = 13, width = 13)
print(ggplot(corPlot, aes(x = Var1, y = Var2, fill = value)) +
    theme_classic() +
    oldham_theme() +
    geom_tile() +
    geom_text(aes(Var1, Var2, label = label), color = "black", size = 7) +
    scale_fill_distiller(palette = "RdBu", breaks = seq(-1, 1, 1), limits = c(-1, 1)) +
    labs(title = "Clonal Correlations", x = "", y = "", fill = "Pearson\ncorrelation") +
    theme(axis.text.x = element_text(angle = 90)))
dev.off()
# create correlation heatmap, noncumulative version
barOutIn = barOut[, -c(1, 2, 4)]
colnames(barOutIn) = c("Clone.1", "Clone.3", "Clone.4", "Clone.5", "Clone.6")
corPlot = reshape2::melt(cor(barOutIn))
corPlot = data.table(corPlot, label = round(corPlot$value, 1))
corPlot$Var1 = gsub("\\.", " ", corPlot$Var1)
corPlot$Var2 = gsub("\\.", " ", corPlot$Var2)
setwd("~/@patrick/SF9495/integration_analysis/glasso/noncumulative")
pdf("cor_heatmap.pdf", height = 13, width = 13)
print(ggplot(corPlot, aes(x = Var1, y = Var2, fill = value)) +
    theme_classic() +
    oldham_theme() +
    geom_tile() +
    geom_text(aes(Var1, Var2, label = label), color = "black", size = 7) +
    scale_fill_distiller(palette = "RdBu", breaks = seq(-1, 1, 1), limits = c(-1, 1)) +
    labs(title = "Noncumulative clonal\nabundance correlations", x = "", y = "", fill = "Pearson\ncorrelation") +
    theme(axis.text.x = element_text(angle = 90)))
dev.off()
