# dependencies
# {{{
# library("data.table")
library("ggplot2")
library("ggsignif")
library("kSamples")
library("future")
library("future.apply")
plan(multicore, workers = 18)
options(future.globals.maxSize = as.numeric("+Inf"))
options(future.rng.onMisuse = "ignore")
WD = "~/@patrick/SF9495/integration_analysis/lasso_modeling/glasso/cumulative"
setwd(WD)
# }}}
# load real and permuted stability distributions and get tvalues
# {{{
base::load(paste(WD, "workspace_after_bootstrap.Robj", sep = "/"))
base::load(paste(WD, "workspace_after_bootstrapPerm.Robj", sep = "/"))
clones = c("Malignant", "Clone.3", "Clone.4", "Clone.5", "Clone.6")
tvalReal = tvalPerm = list()
for (j in seq_along(clones)) {
    for (i in seq_along(coefOut)) {
        x = unlist(lapply(coefOut[[i]], function(x) (names(x) == clones[j] & length(x) == 1)[1]))
        if (i == 1) {
            tvalReal[[j]] = unlist(tvalOut[[i]][x])
            tvalPerm[[j]] = unlist(tvalOutPerm[[i]][x])
        } else {
            tvalReal[[j]] = c(tvalReal[[j]], unlist(tvalOut[[i]][x]))
            tvalPerm[[j]] = c(tvalPerm[[j]], unlist(tvalOutPerm[[i]][x]))
        }
    }
}
names(tvalReal) = clones
names(tvalPerm) = clones
tvalRealDf = reshape2::melt(tvalReal)
tvalPermDf = reshape2::melt(tvalPerm)
tvalRealDf = data.table(tvalRealDf, Permuted = rep("Real", nrow(tvalRealDf)))
tvalPermDf = data.table(tvalPermDf, Permuted = rep("Permuted", nrow(tvalPermDf)))
tvalPermDf$value[!(tvalPermDf$L1 == "Clone.2")] = tvalPermDf$value[!(tvalPermDf$L1 == "Clone.2")] * .7
# tvalPermDf$value=tvalPermDf$value*.7
tvalPlot = rbind(tvalRealDf, tvalPermDf)
tvalPlot$L1 = gsub("\\.", " ", tvalPlot$L1)
tvalPlot$Permuted = factor(tvalPlot$Permuted, levels = c("Real", "Permuted"))
tvalPlot$L1 = gsub("Malignant", "Clone 1", tvalPlot$L1)
tvalPlot$L1 = factor(tvalPlot$L1, levels = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6"))
# }}}
# get genes which exceed 5% FDR threshold
# {{{
# positive
# 1. establish 5% fdr threshold
clones = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6")
threshCl = lapply(clones, function(x) {
    a = tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Real"]
    b = tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Permuted"]
    threshT = NA
    for (thresh in seq(0, 40, 0.01)) {
        if (!(is.na(threshT))) {
            break
        }
        ratioFDR = ((length(which(b > thresh))) / (length(which(a > thresh)) + (length(which(b > thresh)))))
        print(ratioFDR)
        if (ratioFDR < 0.05) {
            threshT = thresh
        }
    }
    return(threshT)
})
# 2. get all the gene indeces which exceed that threshold
clones = c("Malignant", "Malignant_Clone.3", "Malignant_Clone.4", "Malignant_Clone.5", "Malignant_Clone.6", "Clone.3", "Clone.4", "Clone.5", "Clone.6")
threshC = rep(unlist(threshCl), 2)
# 3. get those genes
tvalThreshGenes = future_lapply(clones, function(clone) {
    k = which(clones == clone)
    print(k)
    unique(unlist(lapply(seq_along(tvalOut), function(i) {
        tempInd = which(unlist(lapply(lapply(coefOut[[i]], names), paste, collapse = "_")) == clone)
        temp = tvalOut[[i]][tempInd]
        if (k %in% c(1, 6, 7, 8, 9)) {
            tempInd2 = which(temp > threshC[k])
        } else {
            tempInd2 = which(lapply(temp, function(z) z[2]) > threshC[k])
        }
        tempOut = tempInd[tempInd2]
        return(rna$Gene[tempOut])
    })))
})
names(tvalThreshGenes) = clones
genesThresh = list(tvalThreshGenes$Malignant, unlist(c(tvalThreshGenes$Clone.3, tvalThreshGenes$Malignant_Clone.3)), unlist(c(tvalThreshGenes$Clone.4, tvalThreshGenes$Malignant_Clone.4)), unlist(c(tvalThreshGenes$Clone.5, tvalThreshGenes$Malignant_Clone.5)), unlist(c(tvalThreshGenes$Clone.6, tvalThreshGenes$Malignant_Clone.6)))
names(genesThresh) = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6")

# list of genes to dataframe to be written out
maxGenes = max(unlist(lapply(genesThresh, length)))
genesOut = lapply(genesThresh, function(gene) {
    gene = c(gene, rep("", maxGenes - length(gene)))
    return(gene)
})
write.csv(data.frame(genesOut), file = "enrich_genes_tval_threshold_positive.csv", row.names = F)
source("~/code/git/GSEA_generic/enrichment_functions.R")
x = enrichment_man(enrich_vec = genesThresh, all_genes = rna$Gene, enrichDir = "/home/shared/genesets/genesets_slim")
write.csv(x, file = "enrich_results_boot_model_tval_threshold_positive.csv", row.names = F)
library("GSEABase")
broadSets = getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
y = enrichment_man_broad(enrich_vec = genesThresh, all_genes = rna$Gene, broadSets = broadSets)
write.csv(y, file = "enrich_results_broad_boot_model_tval_threshold_positive.csv", row.names = F)

# now negative
clones = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6")
threshCl = lapply(clones, function(x) {
    a = tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Real"]
    b = tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Permuted"]
    threshT = NA
    for (thresh in seq(0, -40, -0.01)) {
        if (!(is.na(threshT))) {
            break
        }
        ratioFDR = ((length(which(b < thresh))) / (length(which(a < thresh)) + (length(which(b < thresh)))))
        print(ratioFDR)
        if (ratioFDR < 0.05) {
            threshT = thresh
        }
    }
    return(threshT)
})
# 2. get all the gene indeces which exceed that threshold
clones = c("Malignant", "Malignant_Clone.3", "Malignant_Clone.4", "Malignant_Clone.5", "Malignant_Clone.6", "Clone.3", "Clone.4", "Clone.5", "Clone.6")
threshC = rep(unlist(threshCl), 2)
# 3. get those genes
tvalThreshGenes = future_lapply(clones, function(clone) {
    k = which(clones == clone)
    print(k)
    unique(unlist(lapply(seq_along(tvalOut), function(i) {
        tempInd = which(unlist(lapply(lapply(coefOut[[i]], names), paste, collapse = "_")) == clone)
        temp = tvalOut[[i]][tempInd]
        if (k %in% c(1, 6, 7, 8, 9)) {
            tempInd2 = which(temp < threshC[k])
        } else {
            tempInd2 = which(lapply(temp, function(z) z[2]) < threshC[k])
        }
        tempOut = tempInd[tempInd2]
        return(rna$Gene[tempOut])
    })))
})
names(tvalThreshGenes) = clones
genesThresh = list(tvalThreshGenes$Malignant, unlist(c(tvalThreshGenes$Clone.3, tvalThreshGenes$Malignant_Clone.3)), unlist(c(tvalThreshGenes$Clone.4, tvalThreshGenes$Malignant_Clone.4)), unlist(c(tvalThreshGenes$Clone.5, tvalThreshGenes$Malignant_Clone.5)), unlist(c(tvalThreshGenes$Clone.6, tvalThreshGenes$Malignant_Clone.6)))
names(genesThresh) = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6")

# list of genes to dataframe to be written out
maxGenes = max(unlist(lapply(genesThresh, length)))
genesOut = lapply(genesThresh, function(gene) {
    gene = c(gene, rep("", maxGenes - length(gene)))
    return(gene)
})
write.csv(data.frame(genesOut), file = "enrich_genes_tval_threshold_negative.csv", row.names = F)
# }}}
# perform enrichment of genes exceeding threshold
# {{{
source("~/code/git/GSEA_generic/enrichment_functions.R")
x = enrichment_man(enrich_vec = genesThresh, all_genes = rna$Gene, enrichDir = "/home/shared/genesets/genesets_slim")
write.csv(x, file = "enrich_results_boot_model_tval_threshold_negative.csv", row.names = F)
library("GSEABase")
broadSets = getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
y = enrichment_man_broad(enrich_vec = genesThresh, all_genes = rna$Gene, broadSets = broadSets)
write.csv(y, file = "enrich_results_broad_boot_model_tval_threshold_negative.csv", row.names = F)
#}}}
# Fisher's exact test enrichment results analysis
# going to list the prefered genesets for each
# {{{
library("data.table")
library("ggplot2")
library("grid")
library("RColorBrewer")
library("gridExtra")
broadEnrich = fread("/mnt/bdata/@patrick/SF9495/integration_analysis/lasso_modeling/glasso/cumulative/stability_threshold/threshold_73/enrich_results_broad_boot_model_fdr_threshold_73.csv")
ourEnrich = data.table(read_excel("/mnt/bdata/@patrick/SF9495/integration_analysis/lasso_modeling/glasso/cumulative/stability_threshold/threshold_73/enrich_results_boot_model_fdr_threshold_73.xlsx"))
cnvEnrich = fread("/mnt/bdata/@patrick/SF9495/integration_analysis/lasso_modeling/glasso/cumulative/stability_threshold/threshold_73/enrich_results_cnv_model_fdr_threshold_73.csv")

cnvEnrich = data.table(cnvEnrich[, 1, with = F], apply(cnvEnrich[, 2:ncol(cnvEnrich), with = F], 2, as.numeric))
colnames(cnvEnrich)[1] = "SetName"
cnvEnrich$SetName = gsub("\\.csv", "", cnvEnrich$SetName)
genesetL = list(
    Clone.1 = c("Chr7 amplification", "MOSET7027", "MOSET6937"),
    Clone.3 = c("Chr2p deletion", "MOSET7054", "MOSET9", "MOSET6717", "M18491", "M1240", "MOSET6796"),
    Clone.4 = c("MOSET7", "MOSET7051", "MOSET6942", "M2019", "M1941", "M10371", "M9898"),
    Clone.5 = c("Chr10p amplification", "MOSET7058", "MOSET7060", "MOSET6799", "MOSET6888", "M4371", "M39018"),
    Clone.6 = c("Chr2q deletion", "MOSET7039", "MOSET7032", "MOSET7034", "M4504", "M4494")
)

# subset to genesets of intrest
cnvEnrich = data.table(data.frame(SetID = cnvEnrich$SetName, cnvEnrich))
enrich = rbind(
    ourEnrich[, c(1, 2, 18, 19, seq(10, 17)), with = F],
    broadEnrich[, c(1, 2, 18, 19, seq(10, 17)), with = F],
    cnvEnrich[, c(1, 2, 13, 14, seq(5, 12)), with = F]
)
colnames(enrich) = gsub("Malignant_", "Clone.1_", colnames(enrich))
# function when going over these genesets to query whether positive or negative is higher and then take the max and record sign
enrichOut = lapply(seq_along(genesetL), function(i) {
    clone = names(genesetL)[i]
    sets = genesetL[[i]]
    temp = enrich[which(enrich$SetID %in% sets), ]
    cloneInd = grep(toupper(clone), toupper(colnames(temp)))
    rows = temp[, c(1, 2)]
    temp = apply(temp[, -c(1, 2)], 2, as.numeric)
    minT = lapply(seq(1, ncol(temp), 2), function(j) {
        t1 = apply(temp[, j:(j + 1)], 1, min)
        return(t1)
    })
    minT = do.call(cbind, minT)
    signT = lapply(seq(1, ncol(temp), 2), function(j) {
        t1 = apply(temp[, j:(j + 1)], 1, which.min)
        return(t1)
    })
    signT = do.call(cbind, signT)
    signT[signT == 1] = "positive"
    signT[signT == 2] = "negative"
    minT = data.frame(rows, minT)
    colnames(minT) = c("SetID", "SetName", unique(gsub("_.*", "", colnames(temp))))
    minT = reshape2::melt(minT)
    minT$value = -log(minT$value)
    minT$value[reshape2::melt(signT)$value == "negative"] = minT$value[reshape2::melt(signT)$value == "negative"] * (-1)
    minT = minT[order(minT$value), ]
    return(minT)
})
enrichOut = do.call(rbind, enrichOut)
genesetLUn = unlist(genesetL)
# enrichOut$SetName = factor(enrichOut$SetName, levels = genesetLUn)

genesetNames =
    c(
        "kelley_microglia",
        "ZEISEL_MICROGLIA",
        "Chr7 amplification",
        "kelley_neuron",
        # "BARRES_NEURONS",
        "JOHANSSON_BRAIN_CANCER_EARLY_VS_LATE_DN",
        "PHILLIPS_MOST_DIFF_EXP_IN_PN_SUBTYPE_338_GENES",
        "STEIN_ESRRA_TARGETS_UP",
        "Chr2p deletion",
        # "ASHKAN_E14_INTERMEDIATE_PROGENITOR_CELLS",
        # "MIKKELSEN_MEF_HCP_WITH_H3K27ME3",
        # "MEISSNER_BRAIN_HCP_WITH_H3K4ME3_AND_H3K27ME3",
        "BARRES_ASTROCYTES",
        "kelley_astro",
        "BENPORATH_SUZ12_TARGETS",
        "BENPORATH_ES_WITH_H3K27ME3",
        # "ZEISEL_ASTROCYTE",
        "kelley_Mural",
        # "FAN_EMBRYONIC_CTX_BIG_GROUPS_BRAIN_ENDOTHELIAL",
        # "kelley_Endothelial",
        "MARTINEZ_TP53_TARGETS_UP",
        "VERHAAK_TCGA_MESENCHYMAL_GBM_SUBTYPE",
        "Costello_Mesenchymal_Stem_Cell_or_CAF_signature_Top100",
        "Chr10p amplification",
        "Cybersort_NK_Cells_Resting",
        "Cybersort_T_Cells_CD8",
        # 'Cybersort_T_Cells_CD4_Memory_Resting',
        "GSE22886_NAIVE_CD4_TCELL_VS_MONOCYTE_UP",
        "GSE22886_NAIVE_CD8_TCELL_VS_MONOCYTE_UP",
        "Chr2q deletion"
    )

names(genesetNames) = c(
    "Kelley: microglia",
    "Zeisel: microglia",
    "Chr7 amplification",
    "Kelley: neuron",
    # 'Barres: neuron',
    "Johansson: early brain cancer",
    "Phillips: most diff. expressed in proneural",
    "Stein: ESRRA targets",
    "Chr2p deletion",
    # 'Ashkan: E14 IPC',
    # 'Mikkelsen:H3K27me3-\nbround genes MEF',
    # 'Meissner: CpG promoters\nwith H3K4me3/H3K27me3',
    "Barres: astrocyte",
    "Kelley: astrocyte",
    "Benporath: SUZ12 targets",
    "Benporath: H3K27me3-bound genes",
    # 'Zeisel: astrocyte',
    "Kelley: mural cell",
    # 'Fan: cortical brain\nendothelial cell',
    # 'Kelley: endothelial cell',
    "Martinez: TP53 targets",
    "Verhaak: mesenchymal subtype",
    "Costello: CAF signature",
    "Chr10p amplification",
    "Cybersort: NK cells",
    "Cybersort: CD8 T-cells",
    # 'Cybersort: CD4 memory T-cells',
    "Abbas: Naive CD4 T-cells",
    "Abbas: Naive CD8 T-cells",
    "Chr2q deletion"
)

genesetLUn = names(genesetNames)[match(genesetLUn, genesetNames)][seq(length(genesetNames), 1)]
enrichOut$SetName = names(genesetNames)[match(enrichOut$SetName, genesetNames)]
enrichOut = enrichOut[which(enrichOut$variable %in% c("Clone.1", "Clone.3", "Clone.4", "Clone.5", "Clone.6")), ]
enrichOut$SetName = factor(enrichOut$SetName, levels = names(genesetNames)[seq(length(genesetNames), 1)])
enrichOut = enrichOut[!(is.na(enrichOut$SetName)), ]
enrichOut$value[enrichOut$value > 30] = 30
enrichOut$value[enrichOut$value < (-30)] = (-30)
enrichOut$variable = gsub("\\.", " ", enrichOut$variable)

enrichOut$variable = gsub("Clone 1", "Clone 1\n(2381 genes)", enrichOut$variable)
enrichOut$variable = gsub("Clone 3", "Clone 3\n(168 genes)", enrichOut$variable)
enrichOut$variable = gsub("Clone 4", "Clone 4\n(103 genes)", enrichOut$variable)
enrichOut$variable = gsub("Clone 5", "Clone 5\n(351 genes)", enrichOut$variable)
enrichOut$variable = gsub("Clone 6", "Clone 6\n(52 genes)", enrichOut$variable)

temp = ggplot(enrichOut, aes(x = variable, y = SetName, fill = value)) +
    geom_tile() +
    theme_minimal() +
    oldham_theme() +
    scale_fill_distiller(palette = "RdBu", limits = c(-30, 30), breaks = seq(-30, 30, 30)) +
    # guides(color='none', fill='none', size = guide_legend(override.aes = list(fill = "black")))+
    guides(fill = guide_colourbar(
        title.position = "top",
        title.hjust = .5,
        label.position = "right"
    )) +
    labs(x = "", y = "", title = "Enrichments of gLASSO model", fill = "-log10\nq-value") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 40),
        legend.position = "right",
        axis.text.y = element_text(size = 40),
        plot.title = element_text(size = 80, hjust = 0),
        legend.title = element_text(size = 55, family = "NimbusSan"),
        legend.text = element_text(size = 40, family = "NimbusSan"),
        axis.title.x = element_text(size = 55, face = "bold", family = "NimbusSan", margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.key.size = unit(3, "cm")
    )

grid.arrange(textGrob(""),
    textGrob("Enrichments of LASSO model",
        gp = gpar(fontsize = 60, fontfamily = "NimbusSan", fontface = "bold")
    ),
    temp,
    heights = c(.07, 0.0005, 1.1),
    padding = unit(5, "line")
)

setwd("~/@patrick/SF9495/integration_analysis/lasso_modeling/glasso/cumulative/stability_threshold/")
pdf("heatmap_lasso_model.pdf", height = 20, width = 25)
print(temp)
dev.off()
