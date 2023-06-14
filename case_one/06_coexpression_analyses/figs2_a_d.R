# dependencies
# {{{
library("data.table")
library("future")
library("future.apply")
setwd("~/@patrick/SF9495/integration_analysis/network_deconvolution")
library("data.table")
library("forecast")
library("gridExtra")
library("grid")
library("RColorBrewer")
library("ggplot2")
library("egg")
# }}}
# load in clonal frquencies and kME table
# {{{
barOut = fread("~/@patrick/SF9495/figures/sf9495_cellular_abundance.csv")
kME = fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/kME_table_06-26-41.csv")
# }}}
# load RNA
# module we could show here that was significantly enriched with markers of neurons, astrocytes, or mural cells (and not significantly / less significantly associated with clonal abundance compared to the modules featured in Fig. 3)
# {{{
rna = fread("~/@patrick/SF9495/rna_array/Renumbered_SF9495_ALL_69_Qnorm.csv")
meanExpr = future_apply(rna[, -seq(1, 6)], 1, mean)
rna = rna[order(meanExpr, decreasing = T), ]
rna = rna[!(duplicated(rna$Gene)), ]
rna = rna[-grep("^MIR", rna$Gene), ]
rna = rna[-grep("^NCRNA", rna$Gene), ]
rna = rna[-grep("^LOC", rna$Gene), ]
rna = rna[-grep("^SNOR", rna$Gene), ]
rna = rna[-grep("^DKFZ", rna$Gene), ]
rna = rna[-grep("^KIAA", rna$Gene), ]
rna = rna[-grep("^FLJ", rna$Gene), ]
rnaScale = t(future_apply(rna[, -c(seq(1, 6), 69, 70)], 1, log2))
colnames(rnaScale) = gsub("SF9495_", "", colnames(rnaScale))
rnaScale = data.frame(Gene = rna$Gene, rnaScale)
# }}}
# Fisher's exact gene set enrichment
# {{{
# 1 read in enrichments
source("~/code/git/SF10711_clonal_integration_analysis/gseaplot2_oldham.R")
oldhamSets = list(
    fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/oldham_slim_fdr_06-26-41_FDR.csv"),
    fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/oldham_slim_fdr_negative_06-26-41_FDR_NEG.csv")
)
broadSets = list(
    fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/GSHyperG_BROAD_SETS_06-26-41_TOPMODPOSFDR.csv"),
    fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/GSHyperG_BROAD_SETS_06-26-41_FDR_NEG.csv")
)

gtf = fread("/home/shared/hg_align_db/GRCh37.38/gencode.v38lift37.annotation.gene_only.gtf")
gtf = gtf[, -c(2, 3, 6, 8)]
colnames(gtf) = c("chr", "start", "end", "strand", "gene")
gtf$gene = gsub('";.*', "", gsub('.*gene_name "', "", gtf$gene))

genesetL = list(
    green = c("MOSET7053", "MOSET6946", "MOSET6937", "M38974"),
    pink = c("MOSET6955", "MOSET7058", "MOSET6944", "MOSET6956"),
    sienna3 = c("MOSET7052", "MOSET6936", "MOSET8", "MOSET4"),
    greenyellow = c("MOSET9", "MOSET7054", "MOSET7018", "MOSET3")
)

names(genesetL[[2]]) = c("HPA: endothelial\ncells", "Kelley: endothelial\ncells", "LaManno: endothelial\ncells", "GTEX: enodothelial\ncells")
names(genesetL[[3]]) = c("Kelley: oligodendrocytes", "Zeisel: oligodendrocytes", "Barres: oligodendrocytes", "ABA: oligodendrocytes")
names(genesetL[[1]]) = c("Kelley: microglia", "LaManno: microglia", "Zeisel: microglia", "Jinesh: blebbishield\nsignature")
names(genesetL[[4]]) = c("Barres: neurons", "Kelley: neurons", "He: neurons", "ABA: neurons")

# subset to genesets of intrest
oldhamSets = lapply(oldhamSets, function(x) x[which(x$SetID %in% unlist(genesetL)), c(1, 2, which(colnames(x) %in% unlist(c("turquoise", "blue", "black", "midnightblue", names(genesetL))))), with = F])
broadSets = lapply(broadSets, function(x) x[which(x$SetID %in% unlist(genesetL)), c(1, 2, which(colnames(x) %in% unlist(c("turquoise", "blue", "black", "midnightblue", names(genesetL))))), with = F])

# function when going over these genesets to query whether positive or negative is higher and then take the max and record sign
genesetConv = unlist(genesetL)
names(genesetConv) = gsub(".*\\.", "", names(genesetConv))
setEval = function(enrichL, geneS) {
    enrichL = lapply(enrichL, function(x) x[which(x$SetID %in% geneS), ])
    enrichL[[2]] = enrichL[[2]][match(enrichL[[1]]$SetID, enrichL[[2]]$SetID), ]
    enrichL = lapply(enrichL, function(x) melt(x[, -2], id.var = "SetID"))
    enrichL = lapply(enrichL, function(x) {
        x$SetID = names(genesetConv)[match(x$SetID, genesetConv)]
        return(x)
    })
    out = data.frame(matrix(nrow = nrow(enrichL[[1]]), ncol = 2))
    for (i in seq_along(rownames(enrichL[[1]]))) {
        ind = which.min(c(enrichL[[1]]$value[i], enrichL[[2]]$value[i]))
        out[i, 1] = enrichL[[ind]]$value[i]
        if (ind == 1) {
            out[i, 2] = "positive"
        }
        if (ind == 2) {
            out[i, 2] = "negative"
        }
    }
    out = data.frame(enrichL[[1]][, -c(3)], out)
    return(out)
}
# }}}
# create matrix with values and matrix with pos or neg
# {{{
oldhamEval = lapply(genesetL, function(set) setEval(oldhamSets, set))
broadEval = lapply(genesetL, function(set) setEval(broadSets, set))
outEval = lapply(seq_along(oldhamEval), function(i) {
    rbind(oldhamEval[[i]], broadEval[[i]])
})
names(outEval) = names(genesetL)
outEval = lapply(seq_along(outEval), function(i) {
    x = outEval[[i]]
    colnames(x) = c("GeneSet", "Module", "P.value", "Sign")
    x$P.value = as.numeric(x$P.value)
    x$P.value = -log10(x$P.value)
    x$P.value[x$Sign == "negative"] = x$P.value[x$Sign == "negative"] * -1
    # Windsorize values to 30
    x$P.value[x$P.value > 30] = 30
    x$P.value[x$P.value < (-30)] = -30
    x$GeneSet = factor(x$GeneSet, levels = x$GeneSet[x$Module == names(outEval)[i]][order(x$P.value[x$Module == names(outEval)[i]])])
    x$Module = factor(x$Module, levels = c("green", "pink", "sienna3", "greenyellow", "turquoise", "blue", "black", "midnightblue"))
    print(i)
    return(x)
})

mods = c("green", "pink", "sienna3", "greenyellow")
for (i in seq_along(mods)) {
    modsL = unlist(c(mods[i], c("turquoise", "blue", "black", "midnightblue")))
    outEval[[i]] = outEval[[i]][which(outEval[[i]]$Module %in% modsL), ]
}
rnaPlot = data.frame(Gene = rnaScale$Gene, t(apply(rnaScale[, -1], 1, scale)))
# }}}
# create enrichment plot
# {{{
enrichPlotFisher = lapply(outEval, function(x) {
    temp = ggplot(x, aes(x = Module, y = GeneSet, fill = P.value)) +
        geom_tile() +
        theme_minimal() +
        oldham_theme() +
        scale_fill_distiller(palette = "RdBu", limits = c(-30, 30), breaks = seq(-30, 30, 30)) +
        # guides(color='none', fill='none', size = guide_legend(override.aes = list(fill = "black")))+
        guides(fill = guide_colourbar(
            title.position = "top",
            title.hjust = .5,
            label.position = "bottom"
        )) +
        labs(x = "", y = "", title = "", fill = "q-value") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 40),
            legend.position = "bottom",
            axis.text.y = element_text(size = 40),
            plot.title = element_text(size = 60, hjust = 1.2),
            legend.title = element_text(size = 55, family = "NimbusSan"),
            legend.text = element_text(size = 40, family = "NimbusSan"),
            axis.title.x = element_text(size = 55, face = "bold", family = "NimbusSan", margin = margin(t = 0, r = 0, b = 0, l = 0)),
            legend.key.size = unit(3, "cm")
        ) +
        theme(legend.position = "none", plot.margin = margin(b = -35))
})

kmeCols = paste0("kME", mods)
plotLine = lapply(kmeCols, function(kmeCol) {
    plot3 = t(rnaPlot[(rnaPlot$Gene %in% (kME$Gene[order(kME[, which(colnames(kME) == kmeCol), with = F], decreasing = T)[1:12]])), -1])
    colnames(plot3) = rnaScale$Gene[(rnaScale$Gene %in% (kME$Gene[order(kME[, which(colnames(kME) == kmeCol), with = F], decreasing = T)[1:12]]))]
    rownames(plot3) = seq(1, nrow(plot3))
    plot3 = reshape2::melt(plot3)
    p3 = ggplot(plot3, aes(x = Var1, y = value, color = Var2)) +
        geom_line(size = 2) +
        theme_classic() +
        oldham_theme() +
        scale_x_continuous(breaks = c(1, nrow(plot3)), labels = c(1, 81)) +
        scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(12)) +
        guides(color = guide_legend(ncol = 4, byrow = TRUE)) +
        labs(x = "Section ID", y = "Scaled expression", color = "", title = "") +
        theme(
            legend.position = "bottom",
            axis.text.y = element_text(size = 40),
            axis.text.x = element_text(size = 40),
            axis.title.y = element_text(size = 40),
            axis.title.x = element_text(size = 40),
            plot.title = element_text(size = 60),
            legend.text = element_text(size = 32, family = "NimbusSan")
        )
})

pdf("green_mod.pdf", height = 13, width = 30)
grid.arrange(plotLine[[1]], enrichPlotFisher[[1]], ncol = 2, top = textGrob("Green module snapshot", gp = gpar(fontsize = 100, fontfamily = "NimbusSan", fontface = "bold")))
dev.off()

pdf("pink_mod.pdf", height = 13, width = 30)
grid.arrange(plotLine[[2]], enrichPlotFisher[[2]], ncol = 2, top = textGrob("Pink module snapshot", gp = gpar(fontsize = 100, fontfamily = "NimbusSan", fontface = "bold")))
dev.off()

pdf("sienna3_mod.pdf", height = 13, width = 30)
grid.arrange(plotLine[[3]], enrichPlotFisher[[3]], ncol = 2, top = textGrob("Sienna3 module snapshot", gp = gpar(fontsize = 100, fontfamily = "NimbusSan", fontface = "bold")))
dev.off()

pdf("greenyellow_mod.pdf", height = 13, width = 30)
grid.arrange(plotLine[[4]], enrichPlotFisher[[4]], ncol = 2, top = textGrob("Greenyellow module snapshot", gp = gpar(fontsize = 100, fontfamily = "NimbusSan", fontface = "bold")))
dev.off()
# }}}
