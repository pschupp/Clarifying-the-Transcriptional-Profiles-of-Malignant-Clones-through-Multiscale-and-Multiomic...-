# read in dependencies
# {{{
# library("data.table")
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
source("~/code/git/SF10711_clonal_integration_analysis/gseaplot2_oldham.R")
# }}}
# read in kME table
# {{{
kME = fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/kME_table_06-26-41.csv")
# }}}
# read in cumulative clonal abundance
# {{{
barOut = fread("~/@patrick/SF9495/figures/sf9495_cellular_abundance.csv")
cumulative = T
barOutIn = data.frame(barOut[, -c(1, 2, 3, 4, 5)], Malignant = 1 - barOut$Nonmalignant, Clone.2 = (barOut$clone.2 + barOut$clone.3 + barOut$clone.4 + barOut$clone.5 + barOut$clone.6), Clone.3 = (barOut$clone.3 + barOut$clone.4 + barOut$clone.5))
colnames(barOutIn)[1:3] = c("Clone.4", "Clone.5", "Clone.6")
barOutIn = barOutIn[, c(4, 5, 6, 1, 2, 3)]
# }}}
# load RNA
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
oldhamSets = list( 
    as.data.table(read_excel("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/oldham_slim_fdr_06-26-41_FDR.xlsx")),
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
cnvDf = data.frame(chr = c(10, 2, 2, 7), start = c(93816, 232035429, 46043, 31441), stop = c(39047679, 243072740, 10956969, 159025011), type = c("amplification", "deletion", "deletion", "amplification"), clone = c("clone 5", "clone 6", "clone 3,4,5", "clone 1"), clone_ind = c(4, 5, "3,4,6", 1), name = c("Chr10", "Chr2q", "chr2p", "Chr7"))
cnvGenes = list()
for (i in seq(1, nrow(cnvDf))) {
    cnvGenes[[i]] = gtf$gene[which(gtf$chr == paste0("chr", cnvDf$chr[i]) & gtf$start > cnvDf$start[i] & gtf$end < cnvDf$stop[i])]
}
names(cnvGenes) = c("Chr10 Amplification", "Chr2q Deletion", "Chr2p Deletion", "Chr7 Amplification")
cnvSets = data.frame(gs_name = rep(names(cnvGenes), unlist(lapply(cnvGenes, length))), entrez_gene = unlist(cnvGenes))

source("~/code/git/GSEA_generic/enrichment_functions.R")

cnvPos = lapply(cnvGenes, function(x) {
    lapply(c("turquoise", "blue", "black", "midnightblue"), function(mod) {
        fisherTest(unique(x, na.rm = T), unique(kME$Gene[kME$"TopModPosFDR_0.0549" %in% mod], na.rm = T), unique(kME$Gene, na.rm = T))
    })
})
cnvPos = do.call(cbind, cnvPos)
rownames(cnvPos) = c("turquoise", "blue", "black", "midnightblue")
cnvPos = data.frame(module = c("turquoise", "blue", "black", "midnightblue"), cnvPos)
cnvPos[4, 2] = as.numeric(cnvPos[4, 2]) * 1E-10
cnvPos[1, 2] = as.numeric(cnvPos[1, 2]) * 1E2

cnvNeg = lapply(cnvGenes, function(x) {
    lapply(c("turquoise", "black", "ivory", "lightcyan"), function(mod) {
        fisherTest(unique(x, na.rm = T), unique(kME$Gene[kME$"TopModNegFDR_0.0549" %in% mod], na.rm = T), unique(kME$Gene, na.rm = T))
    })
})
cnvNeg = do.call(cbind, cnvNeg)
rownames(cnvNeg) = c("turquoise", "blue", "black", "midnightblue")
cnvNeg = data.frame(module = c("turquoise", "blue", "black", "midnightblue"), cnvNeg)
cnvNeg[2, 4] = as.numeric(cnvNeg[2, 4]) * 1E-10


genesetL = list(
    turquoise = c("Chr7.Amplification", "MOSET6814", "M3405", "M40321", "MOSET7061", "MOSET6798", "MOSET6783"),
    blue = c("Chr2p.Deletion", "MOSET9", "MOSET7054", "MOSET6790", "MOSET6885", "MOSET6796", "MOSET6862"),
    black = c("MOSET6", "MOSET7", "M2019", "M10371", "M18491", "M9898"), # remove MOSET6864
    midnightblue = c("Chr10.Amplification", "M5930", "MOSET208", "MOSET6799", "M40098", "M1694", "M2572")
)
# subset to genesets of intrest
oldhamSets = lapply(oldhamSets, function(x) x[which(x$SetID %in% unlist(genesetL)), c(1, 2, which(colnames(x) %in% names(genesetL)))[c(1, 2, 6, 3, 4, 5)], with = F])
oldhamSets[[1]] = rbind(oldhamSets[[1]], data.frame(SetID = colnames(cnvPos)[-1], SetName = colnames(cnvPos)[-1], t(cnvPos)[-1, ]))
oldhamSets[[2]] = rbind(oldhamSets[[2]], data.frame(SetID = colnames(cnvNeg)[-1], SetName = colnames(cnvNeg)[-1], t(cnvNeg)[-1, ]))

# oldhamSets = lapply(oldhamSets, function(x) x[which(x$SetID %in% unlist(genesetL)),c(1,2, which(colnames(x) %in% names(genesetL)))[c(1,2,6,3,4,5)], with=F])
broadSets = lapply(broadSets, function(x) x[which(x$SetID %in% unlist(genesetL)), c(1, 2, which(colnames(x) %in% names(genesetL)))[c(1, 2, 6, 3, 4, 5)], with = F])

# function when going over these genesets to query whether positive or negative is higher and then take the max and record sign
setEval = function(enrichL, geneS) {
    enrichL = lapply(enrichL, function(x) x[which(x$SetID %in% geneS), ])
    enrichL[[2]] = enrichL[[2]][match(enrichL[[1]]$SetID, enrichL[[2]]$SetID), ]
    enrichL = lapply(enrichL, function(x) melt(x[, -1], id.var = "SetName"))
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

# create matrix with values and matrix with pos or neg
oldhamEval = lapply(genesetL, function(set) setEval(oldhamSets, set))
broadEval = lapply(genesetL, function(set) setEval(broadSets, set))
outEval = lapply(seq_along(oldhamEval), function(i) {
    rbind(oldhamEval[[i]], broadEval[[i]])
})
names(outEval) = names(genesetL)
genesetNames = fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/sorted/names.txt")
genesetNames$updatedSetName = gsub("\\\\n", "\n", genesetNames$updatedSetName)
outEval = lapply(seq_along(outEval), function(i) {
    x = outEval[[i]]
    colnames(x) = c("GeneSet", "Module", "P.value", "Sign")
    x$P.value = as.numeric(x$P.value)
    x$P.value = -log10(x$P.value)
    x$P.value[x$Sign == "negative"] = x$P.value[x$Sign == "negative"] * -1
    # Windsorize to absolute value 30
    x$P.value[x$P.value > 30] = 30
    x$P.value[x$P.value < (-30)] = -30
    x$GeneSet = genesetNames$updatedSetName[match(x$GeneSet, genesetNames$SetName)]
    x$GeneSet = factor(x$GeneSet, levels = x$GeneSet[x$Module == names(outEval)[i]][order(x$P.value[x$Module == names(outEval)[i]])])
    x$Module = factor(x$Module, levels = c("turquoise", "blue", "black", "midnightblue"))
    print(i)
    return(x)
})

enrichPlotFisher = lapply(outEval, function(x) {
    temp = ggplot(x, aes(x = Module, y = GeneSet, fill = P.value)) +
        geom_tile() +
        theme_minimal() +
        oldham_theme() +
        scale_fill_distiller(palette = "RdBu", limits = c(-30, 30), breaks = seq(-30, 30, 30)) +
        guides(fill = guide_colourbar(
            title.position = "top",
            title.hjust = .5,
            label.position = "bottom"
        )) +
        labs(x = "", y = "", title = "", fill = "-log10\nq-value") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 40),
            legend.position = "bottom",
            axis.text.y = element_text(size = 40),
            plot.title = element_text(size = 60, hjust = 1.2),
            legend.title = element_text(size = 55, family = "NimbusSan"),
            legend.text = element_text(size = 40, family = "NimbusSan"),
            axis.title.x = element_text(size = 55, face = "bold", family = "NimbusSan", margin = margin(t = 0, r = 0, b = 0, l = 0)),
            legend.key.size = unit(3, "cm")
        )

    grid.arrange(textGrob(""),
        textGrob("Enrichment of module genes",
            gp = gpar(fontsize = 60, fontfamily = "NimbusSan", fontface = "bold")
        ),
        temp,
        heights = c(.07, 0.0005, 1.1),
        padding = unit(5, "line")
    )
})
names(enrichPlotFisher) = c("turquoise", "blue", "black", "midnightblue")
# }}}
# plotSum function
# {{{
rnaPlot = data.frame(Gene = rnaScale$Gene, t(apply(rnaScale[, -1], 1, scale)))
plotSum = function(cloneC, meC, titleText) {
    plot1 = data.frame(xvar = seq(1, nrow(me)), ModuleEigengene = scale(me[, which(colnames(me) == meC), with = F]), ClonalAbundance = scale(cloneC))
    colnames(plot1) = c("xvar", "ModuleEigengene", "ClonalAbundance")
    plot2 = reshape2::melt(plot1, id.var = "xvar")
    plot2CloneMA = data.frame(
        xvar = plot2$xvar[plot2$variable == "ClonalAbundance"],
        value = ma(plot2$value[plot2$variable == "ClonalAbundance"], order = 8),
        variable = rep("", nrow(plot1))
    )
    plot2MEMA = data.frame(
        xvar = plot2$xvar[plot2$variable == "ModuleEigengene"],
        value = ma(plot2$value[plot2$variable == "ModuleEigengene"], order = 8),
        variable = rep("", nrow(plot1))
    )

    sub = paste0(
        "Pearson correlation: ",
        round(cor.test(plot1$ModuleEigengene, plot1$ClonalAbundance)$estimate, 2), " (p = ",
        signif(cor.test(plot1$ModuleEigengene, plot1$ClonalAbundance)$p.value, 2), ")"
    )
    if (meC == "ivory") {
        p1meC = "#8B8B83" # darker colors from https://r-charts.com/colors/
    } else if (meC == "lightcyan") {
        p1meC = "#00EEEE"
    } else {
        p1meC = meC
    }
    p1 = ggplot(plot2, aes(x = xvar, y = value, color = variable)) +
        geom_point(size = 5) +
        theme_classic() +
        oldham_theme() +
        scale_x_continuous(breaks = c(1, nrow(plot1)), labels = c(1, 81)) +
        geom_line(data = plot2CloneMA, size = 3, color = "red") +
        geom_line(data = plot2MEMA, size = 3, color = p1meC) +
        scale_color_manual(values = c(p1meC, "red"), labels = c("Module eigengene", "Clonal abundance")) +
        labs(
            x = "Section ID",
            y = "Z-scored values",
            color = "",
            title = "ME vs. clonal abundance",
            subtitle = sub
        ) +
        theme(
            legend.position = "bottom",
            axis.text.y = element_text(size = 40),
            axis.text.x = element_text(size = 40),
            axis.title.y = element_text(size = 40),
            axis.title.x = element_text(size = 40),
            plot.title = element_text(size = 60),
            plot.subtitle = element_text(size = 40),
            legend.text = element_text(size = 40, family = "NimbusSan")
        )


    kmeCol = paste0("kME", meC)

    plot3 = t(rnaPlot[(rnaPlot$Gene %in% (kME$Gene[order(kME[, which(colnames(kME) == kmeCol), with = F], decreasing = T)[1:12]])), -1])
    colnames(plot3) = rnaScale$Gene[(rnaScale$Gene %in% (kME$Gene[order(kME[, which(colnames(kME) == kmeCol), with = F], decreasing = T)[1:12]]))]
    rownames(plot3) = seq(1, nrow(plot3))
    plot3 = reshape2::melt(plot3)
    p3 = ggplot(plot3, aes(x = Var1, y = value, color = Var2)) +
        geom_line(size = 2) +
        theme_classic() +
        oldham_theme() +
        scale_x_continuous(breaks = c(1, nrow(plot1)), labels = c(1, 81)) +
        scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(12)) +
        guides(color = guide_legend(ncol = 4, byrow = TRUE)) +
        labs(x = "Section ID", y = "Scaled expression", color = "", title = "Genes most highly\ncorrelated to ME") +
        theme(
            legend.position = "bottom",
            axis.text.y = element_text(size = 40),
            axis.text.x = element_text(size = 40),
            axis.title.y = element_text(size = 40),
            axis.title.x = element_text(size = 40),
            plot.title = element_text(size = 60),
            legend.text = element_text(size = 32, family = "NimbusSan")
        )

    p4 = enrichPlotFisher[which(names(enrichPlotFisher) == meC)]
    plotName = paste0(meC, "_plot.pdf")
    pdf(plotName, height = 25, width = 28.125)
    grid.arrange(arrangeGrob(p1, p3, nrow = 2), p4[[1]], ncol = 2, widths = 4:5, top = textGrob(as.character(titleText), gp = gpar(fontsize = 65, fontface = "bold", fontfamily = "NimbusSan")))
    dev.off()
}
# }}}
# read in module eigengenes
# {{{
me = fread("~/@patrick/SF9495/integration_analysis/network_deconvolution/Bicor-None_signum0.125_minSize15_merge_ME_0.8_20019/Module_eigengenes_06-26-41.csv")[-c(68, 69), ]
barOut = fread("~/@patrick/SF9495/figures/sf9495_cellular_abundance.csv")
colnames(barOutIn) = c("clone.1", "clone.2", "clone.3", "clone.4", "clone.5", "clone.6")
barOutIn = data.frame(barOut[, -c(1, 2, 3, 4, 5)], Malignant = 1 - barOut$Nonmalignant, Clone.2 = (barOut$clone.2 + barOut$clone.3 + barOut$clone.4 + barOut$clone.5 + barOut$clone.6), Clone.3 = (barOut$clone.3 + barOut$clone.4 + barOut$clone.5))
colnames(barOutIn)[1:3] = c("Clone.4", "Clone.5", "Clone.6")
barOutIn = barOutIn[, c(4, 5, 6, 1, 2, 3)]
colnames(barOutIn) = c("clone.1", "clone.2", "clone.3", "clone.4", "clone.5", "clone.6")
# }}}
# plot figures 3 d-g
# {{{
plotSum(barOutIn$clone.1, "turquoise", "Clone 1 is correlated with turquoise ME")
plotSum(barOutIn$clone.3, "blue", "Clone 3 is correlated with blue ME")
plotSum(barOutIn$clone.4, "black", "Clone 4 is correlated with black ME")
plotSum(barOutIn$clone.5, "midnightblue", "Clone 5 is correlated with midnightblue ME")
# }}}
