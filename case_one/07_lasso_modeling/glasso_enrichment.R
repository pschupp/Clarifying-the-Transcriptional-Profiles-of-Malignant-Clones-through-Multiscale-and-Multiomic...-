
# first look ath the regular and permuted coefficients and make panel for each clone. Use violin plot. Get tvalues...don't not get these...will have to redo or convert
# goal of the scriptis to do enrichment of lasso modelling
setwd("~/@patrick/SF9495/integration_analysis/glasso")
base::load("~/@patrick/SF9495/integration_analysis/glasso/workspace_after_bootstrap.Robj")
base::load("~/@patrick/SF9495/integration_analysis/glasso/workspace_after_bootstrapPerm.Robj")

# to prevent confusion will label t he following variables which have been imported in two flavors, regular and from permuted data
# coefOut - length 100, for all bootstraps coefficients for model with minimum AIC, this model is the relaxed lasso model, so no restriction on scale of coefficients
# pvalOut - length 100, for all bootstraps model p-value for model with minimum AIC
# tvalsOut - for all bootstraps, t-values for all coefficients in the model, length will the same for each member as coefOut
# lengthBag - average number of coefficients for each gene across all bootstraps
# stability - highest indcidence of a combination of clone vector (independent variable) across all bootsraps
# pvalBag - p-value across all bootstraps for each gene averaged using Fisher's method

# look at 95/99th precentile of stability for real v. permuted and use those models only and then go gene enrichment
# also make figure of t-values

data.frame(real = quantile(stability, seq(0, 1, .05), na.rm = T), permuted = quantile(stabilityPerm, seq(0, 1, .05), na.rm = T))

quantile(stabilityPerm, seq(.9, 1, .01))

data.frame(real = quantile(pvalBag, seq(0, 1, .05), na.rm = T), permuted = quantile(pvalBagPerm, seq(0, 1, .05), na.rm = T))
library("future")
library("future.apply")
library("ggplot2")
library("qvalue")
plan(multicore, workers = 18)
options(future.globals.maxSize = as.numeric("+Inf"))
options(future.rng.onMisuse = "ignore")

# mode (most prevalent) factor combination for the model
bootMode = future_lapply(seq(1, nrow(rnaScale)), function(i) names(which.max(table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse = "_")))))))
x = sort(table(unlist(bootMode[stability > 72])), decreasing = T)
# clone 2 is fundemtnatlly invariant. hard to find genes that follow this expression pattern. Clone, suprinsginlgy we can find many patterns, is this model using the cumulative or noncumulative clonal definitions? We ARE using the cumulative defintions. even then we are findding manuy genes uniuqely best modeled by clone 1. Suggests this approach might be less reliant on variantion than FM, though both are largeley based on variance/covariance. But FM has difficult adjusting for low variacne against a background of high variance, i.e. it is optmized to find the most salient features, where as a correlative model, especially one using bootstrapping like this might be able to find the signal despite the large dynamic range of variations.
# next step, knowing the bootMode, use that to subset bagged pval calculations to only those using the mode. after that try replotting the graph and make a cutoff, presumably around 28-35 in stability
# first, for each gene, select only the bootstraps which are the boot mode
modeSelec = future_lapply(seq(1, nrow(rnaScale)), function(i) which(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse = "_"))) == bootMode[i]))
# now need to use modeSelec to only sample from the correct bootstraps
future_lapply(seq(1, nrow(rnaScale)), function(i) x[[1]])
pvalBagMode = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) mean(unlist(lapply(pvalOut[modeSelec[[i]]], function(x) x[[i]])))))
pvalBagNoMode = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) mean(unlist(lapply(pvalOut[-modeSelec[[i]]], function(x) x[[i]])))))

plotDf = data.frame(as.factor(stability), pvalBagMode, pvalBagNoMode, lengthBag)
p = list()
p[[1]] = ggplot(plotDf, aes(x = as.factor(stability), y = (log(pvalBagMode + 1E-17) * -1))) +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme() +
    scale_x_discrete(breaks = seq(0, 100, 20)) +
    labs(x = "Stability", y = "-log10(p-value)", title = "Mode bootstrap: p-value\nincreases as stability increases")
p[[2]] = ggplot(plotDf, aes(x = as.factor(stability), y = (log(pvalBagNoMode + 1E-17) * -1))) +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme() +
    scale_x_discrete(breaks = seq(0, 100, 20)) +
    labs(x = "Stability", y = "-log10(p-value)", title = "Nonmode bootstrap: p-value\nincreases as stability increases")
p[[3]] = ggplot(plotDf, aes(x = as.factor(round(lengthBag)), y = (log(pvalBagMode + 1E-17) * -1))) +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme() +
    labs(x = "Number of factors", y = "-log10(p-value)", title = "Mode bootstrap: p-value\nas number of factors increases")
p[[4]] = ggplot(plotDf, aes(x = as.factor(round(lengthBag)), y = (log(pvalBagNoMode + 1E-17) * -1))) +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme() +
    labs(x = "Number of factors", y = "-log10(p-value)", title = "Nonmode bootstrap: p-value\nas number of factors increases")

pdf("bagged_values_mode_cumulative.pdf", height = 13, width = 13)
print(p)
dev.off()

# stability greater than 72
# use fdr-pvalue
# get genes with +/- correlations, do enrichment
stabThresh = 72 
x = sort(table(unlist(bootMode[stability > stabThresh & qvalue(pvalBag)$qvalue < 0.05])), decreasing = T)
# get list of genes with positive and negative coefficients
bootGenes = list()
i = 1
j = 1
cloneVec = sort(unique(unlist(strsplit(unique(unlist(bootMode)), "_"))))
for (cloneName in cloneVec) {
    temp1 = which(bootMode[stability > stabThresh & qvalue(pvalBag)$qvalue < 0.05] == cloneName | bootMode[stability > stabThresh & qvalue(pvalBag)$qvalue < 0.05] == paste0("Malignant_", cloneName))
    temp2 = seq(1, nrow(rnaScale))[stability > stabThresh & qvalue(pvalBag)$qvalue < 0.05][temp1]
    temp3 = apply(rnaScale[temp2, ], 1, function(rnaRow) lm(rnaRow ~ barOutIn[, which(colnames(barOutIn) == cloneName)]))
    temp4 = unlist(lapply(temp3, function(x) x$coefficients[2]))
    t4pos = which(temp4 > 0)
    t4neg = which(temp4 < 0)
    genesP = rna$Gene[temp2][t4pos]
    genesN = rna$Gene[temp2][t4neg]
    bootGenes[[j]] = genesP
    bootGenes[[j + 1]] = genesN
    j = j + 2
    i = i + 1
}
lapply(bootGenes, length)
# fisher's exact test enrichment calculation of oldham and broad sets
# {{{
names(bootGenes) = paste(rep(cloneVec, each = 2), rep(c("_positive", "_negative"), times = length(cloneVec)), sep = "")

source("~/code/git/GSEA_generic/enrichment_functions.R")
x = enrichment_man(enrich_vec = bootGenes, all_genes = rna$Gene, enrichDir = "/home/shared/genesets/genesets_slim")
write.csv(x, file = "enrich_results_boot_model_fdr_threshold.csv", row.names = F)
library("GSEABase")
broadSets = getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
y = enrichment_man_broad(bootGenes, rna$Gene, broadSets = broadSets)
write.csv(y, file = "enrich_results_broad_boot_model_fdr_threshold.csv", row.names = F)
# }}}
# Fisher's exact test enrichment result of CNV sets
# {{{
library("data.table")
gtf = fread("/home/shared/hg_align_db/GRCh37.38/gencode.v38lift37.annotation.gene_only.gtf")
gtf = gtf[, -c(2, 3, 6, 8)]
colnames(gtf) = c("chr", "start", "end", "strand", "gene")
gtf$gene = gsub('";.*', "", gsub('.*gene_name "', "", gtf$gene))

cnvDf = data.frame(
    chr = c(2, 18, 11, 7, 8, 9),
    start = c(98204464, 41955206, 55838392, 192969, 11977, 12705),
    stop = c(224844864, 77927112, 134605703, 159144957, 146289882, 141124247),
    type = c("deletion", "deletion", "amplification", "amplification", "amplification", "amplification"),
    clone = c("Clone 3", "Clone 3", "Clone 5", "Clone 1", "Clone 1", "Clone 1"),
    clone_ind = c(3, 3, 5, 1, 1, 1),
    names = c("Chr2q", "Chr18", "Chr11q", "Chr7", "Chr8", "Chr9")
)
cnvGenes = list()
for (i in seq(1, nrow(cnvDf))) {
    cnvGenes[[i]] = gtf$gene[which(gtf$chr == paste0("chr", cnvDf$chr[i]) & gtf$start > cnvDf$start[i] & gtf$end < cnvDf$stop[i])]
}
names(cnvGenes) = c("Chr2q deletion", "Chr18 deletion", "Chr11q amplification", "Chr7 amplification", "Chr8 amplification", "Chr9 amplification")
cnvSets = data.frame(gs_name = rep(names(cnvGenes), unlist(lapply(cnvGenes, length))), entrez_gene = unlist(cnvGenes))
# get the cnvGS genes into files so we can run enrichment_man
setwd("~/@patrick/SF10711/figures/fig6/glasso/cnv_gene_set")
for (i in seq_along(cnvGenes)) {
    write.table(cnvGenes[[i]], paste0(names(cnvGenes)[i], ".csv"), row.names = F, quote = F, col.names = F, sep = ",")
}
names(bootGenes) = paste(rep(cloneVec, each = 2), rep(c("_positive", "_negative"), times = length(cloneVec)), sep = "")
source("~/code/git/GSEA_generic/enrichment_functions.R")
cnvEnrich = enrichment_man(bootGenes, rna$Gene, "~/@patrick/SF10711/figures/fig6/glasso/cnv_gene_set")
setwd("~/@patrick/SF10711/figures/fig6/glasso")
write.table(cnvEnrich, file = "cnvEnrichment_results_boot.csv", row.names = F, sep = ",", quote = F)
# GSEA enrichment calculation of oldham and broad sets
# {{{
library("clusterProfiler")
library("enrichplot")
WD = getwd()
setwd("~/code/GSEA/genesets_slim")
myLegend = read.csv(list.files()[grep("^MyGeneSetsLEGEND\\.csv$", list.files())])
mySetNames = list.files()[grep("MOSET", list.files())]
ordervec = gsub("MOSET", "", mySetNames)
ordervec = gsub(".csv", "", ordervec)
mySetNames = mySetNames[order(as.numeric(ordervec))]
mySets = vector(mode = "list", length = length(mySetNames))
for (i in c(1:length(mySetNames))) {
    mySets[[i]] = read.csv(mySetNames[i], header = F)
}
mySets = sapply(mySets, unlist)
mySets = sapply(mySets, as.character)
mySets = sapply(mySets, toupper)
names(mySets) = myLegend$SetName
mySets = data.frame(gs_name = rep(names(mySets), unlist(lapply(mySets, length))), entrez_gene = unlist(mySets))
setwd(WD)

bootGenes = list()
i = 1
j = 4
bootGseaOut = data.frame(matrix(nrow = length(unique(mySets$gs_name)), ncol = 13))
cloneVec = sort(unique(unlist(strsplit(unique(unlist(bootMode)), "_"))))
for (cloneName in cloneVec) {
    temp1 = which(bootMode[stability > 72] == cloneName)
    temp2 = seq(1, nrow(rnaScale))[stability > 72][temp1]
    temp3 = apply(rnaScale[temp2, ], 1, function(rnaRow) lm(rnaRow ~ barOutIn$Malignant + barOutIn))
    temp4 = unlist(lapply(temp3, function(x) x$coefficients[3]))
}
names(temp4) = c(rna$Gene[temp2])
add = sample(seq(-1E-5, 1E-5, length = nrow(rnaScale)), nrow(rnaScale) - length(temp2))
#    add=rep(0, (nrow(rnaScale)-length(temp2)))
names(add) = rna$Gene[-temp2]
temp4 = c(temp4, add)
#    temp4=unlist(scale(temp4)[,1])
enrich = GSEA(geneList = sort(temp4, decreasing = T), exponent = 1, TERM2GENE = mySets, by = "fgsea", verbose = T, pvalueCutoff = 1, eps = 1e-50, scoreType = "std", pAdjustMethod = "fdr", maxGSSize = 2000)
results = enrich@result
rownames(results) = seq(1, nrow(results))
head(results[, -c(10, 11, 12)], 20)
if (j == 1) {
    bootGseaOut[, 1:3] = data.frame(results$ID, results$NES, results$pvalue)
} else {
    results = results[match(bootGseaOut[, 1], results$ID), ]
    bootGseaOut[, j:(j + 1)] = data.frame(results$NES, results$pvalue)
}
# }}}
