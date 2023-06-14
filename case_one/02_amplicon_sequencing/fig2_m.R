# generating stacked barplot of clonal abundances
# read in data and tidy up
# {{{
library("stringr")
oldham_theme = function() {
    theme(
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 40, color = "black", family = "NimbusSan", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "NimbusSan", margin = margin(r = 10)),
        axis.title.y = element_text(size = 40, face = "bold", family = "NimbusSan", margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.title.x = element_text(size = 40, face = "bold", family = "NimbusSan", margin = margin(t = -10, r = 0, b = 0, l = 0)),
        plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(t = -20, b = 0)),
        plot.subtitle = element_text(size = 40, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(t = 40, b = -40)),
        axis.line.x = element_line(size = 3),
        axis.line.y = element_line(size = 3),
        plot.margin = unit(c(4, 2, 1, 2), "lines"),
        legend.position = "right",
        legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size = 30, family = "NimbusSan")
    )
}
setwd("~/@patrick/SF9495/figures")
cluster = data.frame(fread("~/@patrick/SF9495/clonality/pyclone/outputs/tables/cluster.tsv"))
cluster = reshape(cluster[, c(1, 2, 4)], idvar = "sample_id", timevar = "cluster_id", direction = "wide")
rownames(cluster) = cluster[, 1]
cluster = cluster[, -1]
cluster[is.na(cluster)] = .99
amp = read.table("~/@patrick/SF9495/ampseq/Combined_variant_frequencies_run1_and_run2_with_readCounts.csv", sep = ",", header = T)
amp = amp[-(nrow(amp)), ]
amp$sample = paste("sf", str_pad(gsub(".*_", "", amp$sample), width = 3, pad = 0, side = "left"), sep = "")
cluster = cluster[order(rownames(cluster)), ]
amp = amp[match(rownames(cluster), amp$sample), ]
barOut = fread("~/@patrick/SF9495/clonality/pyclone_cITUP/outputs/tmp/tmp/tree/7/results", skip = "#clone_freq")
colnames(barOut) = c("Clone 1", "Clone 2", "Clone 3", "Clone 4", "Clone 5", "Clone 6", "Clone 7")
barOut = data.frame(section = seq(1, nrow(barOut)), barOut[, -1])
barOut = barOut[, c(1, 7, 2, 4, 5, 6, 3)]
colnames(barOut) = c("section", "clone 1", "clone 2", "clone 3", "clone 4", "clone 5", "clone 6")
# }}}
# corrective scaling of cITUP output to account for tumor purity
# {{{
scaleFac = (amp$IDH1.242.GtoA_freq * 2) / apply(barOut[, -1], 1, sum)
cluster = data.frame(apply(barOut[, -1], 2, function(x) (x * scaleFac)))
barOut = data.frame(section = barOut[, 1], Nonmalignant = 1 - (amp$IDH1.242.GtoA_freq * 2), cluster)
barOutP = melt(data.frame(section = barOut[, 1], Nonmalignant = 1 - (amp$IDH1.242.GtoA_freq * 2), cluster), id.var = "section")
# }}}
# create stacked barplot
# {{{
barOutP$section = factor(barOutP$section, levels = unique(barOutP$section))
barOutP$section = as.numeric(barOutP$section)
pdf("clonal_abundance.pdf", width = 24, height = 13)
print(ggplot(barOutP, aes(x = section, y = value, color = variable, fill = variable)) +
    theme_classic() +
    oldham_theme() +
    geom_bar(stat = "identity", width = 1, size = 2) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 1), labels = c(0, 100)) +
    scale_x_continuous(expand = c(0, 0), limits = c(1, 67), breaks = c(1, 67), labels = c(1, 81)) +
    labs(fill = "Clone", title = "Cellular fraction", x = "Section ID", y = "Percentage") +
    scale_fill_manual(values = c("white", "#ff7f00", "#ffbc79", "#dd95e8", "#984ea3", "#4daf40", "#b15928"), name = "Clones", labels = c("Non-malignant", "Clone 1", "Clone 2", "Clone 3", "Clone 4", "Clone 5", "Clone 6")) +
    scale_color_manual(values = c("black", "#ff7f00", "#ffbc79", "#dd95e8", "#984ea3", "#4daf40", "#b15928"), name = "Clones", labels = c("Non-malignant", "Clone 1", "Clone 2", "Clone 3", "Clone 4", "Clone 5", "Clone 6")))
dev.off()

kme = data.frame(rna[, seq(1, 6), with = F], cor(t(rna[, -c(seq(1, 6), 68, 69)]), barOut[, -1]))
colnames(kme)[seq(7, 13)] = c("kME_Parental", "kME_Clone1", "kME_Clone2", "kME_Clone3", "kME_Clone4", "kME_Clone5", "kME_Clone6")
write.csv(kme, file = "sf9495_kme_table.csv", row.names = F)
write.csv(barOut, file = "sf9495_cellular_abundance.csv", row.names = F)
# }}}
