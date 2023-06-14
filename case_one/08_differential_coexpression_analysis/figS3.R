# moving forward with edgington p-value combination. qvalue results are similar to before.
# 1) subset both expression matrices to the shared genes
# {{{
sub9495 = rnaScale9495[rnaScale9495$Gene %in% geneL$sigPosFDR, ]
sub10711 = rnaScale10711[rnaScale10711$Gene %in% geneL$sigPosFDR, ]
# }}}
# 2) regress out malignant vector for each
# {{{
resid9495 = apply(sub9495[, -1], 1, function(x) residuals(lm(as.numeric(x) ~ barOutIn9495$Malignant)))
resid10711 = apply(sub10711[, -1], 1, function(x) residuals(lm(as.numeric(x) ~ barOutIn10711$Malignant)))
# }}}
# 3) compute correlation vector for each
# {{{
corResid9495 = cor(resid9495)
corResid10711 = cor(resid10711)
# }}}
# 4) take average cor at each location
# {{{
meanCor = median(as.matrix(corResid9495[1:5, 1:5]), as.matrix(corResid10711[1:5, 1:5]))
meanCor = list(corResid9495, corResid10711)
meanCorA = do.call(cbind, meanCor)
meanCorA = array(meanCorA, dim = c(dim(meanCor[[1]]), length(meanCor)))
meanCorOut = future_apply(meanCorA, c(1, 2), mean)
# }}}
# 5) perform FM on megaset with average correlation matrix, color by each turn/sample change
# {{{
source("~/code/git/FindModules/FindModules.lint.par.R")
# source('~/code/git/FindModules/FindModules.R')
expr = data.frame(Gene = geneL$sigPosFDR, t(resid10711))
meanCorOut[diag(meanCorOut)] = 0
meanCorOut = meanCorOut[order(expr$Gene), order(expr$Gene)]
rownames(meanCorOut) = colnames(meanCorOut) = expr$Gene
expr = expr[order(expr$Gene), ]
setwd("/home/patrick/code/git/integration_SF10711_SF9495")
FindModules(
    projectname = "averageCorlow",
    expr = expr,
    geneinfo = 1,
    sampleindex = seq(2, ncol(expr)),
    samplegroups = as.factor(c(rep("A", 40), rep("B", 45))),
    subset = NULL,
    simMat = meanCorOut,
    saveSimMat = FALSE,
    simType = "Pearson",
    overlapType = "None",
    TOtype = "signed",
    TOdenom = "min",
    beta = 1,
    MIestimator = "mi.mm",
    MIdisc = "equalfreq",
    signumType = "rel",
    iterate = TRUE,
    signumvec = seq(0.6, 0.8, 0.1),
    minsizevec = c(8, 10, 12, 15, 20),
    signum = NULL,
    minSize = NULL,
    merge.by = "ME",
    merge.param = 0.85,
    export.merge.comp = FALSE,
    ZNCcut = 2,
    calcSW = FALSE,
    loadTree = FALSE,
    writeKME = TRUE,
    calcBigModStat = FALSE,
    writeModSnap = TRUE
)
# }}}
# 6) do enrichment for GSEA as well as stringdb
# try with string db biggest module
# {{{
mod = read.csv("~/code/git/integration_SF10711_SF9495/averageCorlow_Modules/Pearson-None_signum0.219_minSize8_merge_ME_0.85_3069/kME_table_08-44-38.csv")
modGenes = lapply(unique(mod$TopModPosFDR_0.257), function(x) mod$Gene[which(mod$TopModPosFDR_0.257 == x)])
modGenes = modGenes[!is.na(unique(mod$TopModPosFDR_0.257))]
pdf("mod_string_hits.pdf")
lapply(modGenes, function(x) {
    example1_mapped <- string_db$map(data.frame(gene = x), "gene", removeUnmappedRows = TRUE)
    options(SweaveHooks = list(fig = function() {
        par(mar = c(2.1, 0.1, 4.1, 2.1))
    }))
    hits <- example1_mapped$STRING_id
    getOption("SweaveHooks")[["fig"]]()
    string_db$plot_network(hits)
})
dev.off()
modGenes = lapply(unique(mod$TopModPosBC_1.36e.06), function(x) mod$Gene[which(mod$TopModPosBC_1.36e.06 == x)])
modGenes = modGenes[!is.na(unique(mod$TopModPosBC_1.36e.06))]
pdf("mod_bc_string_hits.pdf")
lapply(modGenes, function(x) {
    example1_mapped <- string_db$map(data.frame(gene = x), "gene", removeUnmappedRows = TRUE)
    options(SweaveHooks = list(fig = function() {
        par(mar = c(2.1, 0.1, 4.1, 2.1))
    }))
    hits <- example1_mapped$STRING_id
    getOption("SweaveHooks")[["fig"]]()
    string_db$plot_network(hits)
})
dev.off()
modGenes = lapply(unique(mod$ModSeed), function(x) mod$Gene[which(mod$ModSeed == x)])
modGenes = modGenes[!is.na(unique(mod$ModSeed))]
pdf("mod_seed_string_hits.pdf")
lapply(modGenes, function(x) {
    example1_mapped <- string_db$map(data.frame(gene = x), "gene", removeUnmappedRows = TRUE)
    options(SweaveHooks = list(fig = function() {
        par(mar = c(2.1, 0.1, 4.1, 2.1))
    }))
    hits <- example1_mapped$STRING_id
    getOption("SweaveHooks")[["fig"]]()
    string_db$plot_network(hits)
})
dev.off()
# }}}
