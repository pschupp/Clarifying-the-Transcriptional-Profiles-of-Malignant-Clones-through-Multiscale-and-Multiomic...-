# violin plot figure
# dependencies
# {{{
library("data.table")
library("ggplot2")
library("ggsignif")
library("kSamples")
# }}}
# normal lasso
# {{{
WD = "~/@patrick/SF9495/integration_analysis/lasso_modeling/glasso/noncumulative"
setwd(WD)
base::load(paste(WD, "workspace_after_bootstrap.Robj", sep = "/"), verbose = TRUE)
base::load(paste(WD, "workspace_after_bootstrapPerm.Robj", sep = "/"), verbose = TRUE)
clones = c("Clone.1", "Clone.3", "Clone.4", "Clone.5", "Clone.6")
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
tvalPlot = rbind(tvalRealDf, tvalPermDf)
tvalPlot$L1 = gsub("\\.", " ", tvalPlot$L1)
tvalPlot$Permuted = factor(tvalPlot$Permuted, levels = c("Real", "Permuted"))
tvalPlot$L1 = gsub("Clone.1", "Clone 1", tvalPlot$L1)
tvalPlot$L1 = factor(tvalPlot$L1, levels = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6"))
clones = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6")
pval = lapply(clones, function(x) {
    print(x)
    samp3k = seq(1, length(tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Real"]))
    z = ad.test(tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Real"][samp3k], tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Permuted"][samp3k])$ad[2, 3]
    return(z)
})
p = ggplot(tvalPlot, aes(x = L1, y = value, fill = Permuted)) +
    geom_violin() +
    theme_classic() +
    oldham_theme() +
    scale_y_continuous(limits = c(-45, 45), breaks = seq(-45, 45, 15)) +
    labs(x = "", y = "T-values", fill = "", title = "T-values of lasso model") + # title='T-values of cumulative lasso model', subtitle='P-values from two-sided Anderson-Darling test') +
    scale_fill_manual(values = c("black", "white")) +
    theme(
        legend.spacing.y = unit(1, "cm"),
        axis.text.x = element_text(angle = 90)
    ) +
    geom_signif(
        y_position = c(rep(35, 5)),
        xmin = seq(0.8, 4.8, 1),
        xmax = seq(1.2, 5.2, 1),
        annotation = signif(unlist(pval), 1),
        tip_length = 0,
        size = 2,
        textsize = 10,
        vjust = -.6,
        family = "NimbusSan"
    ) +
    guides(fill = guide_legend(byrow = TRUE))

pdf("tvalue_violin.pdf", height = 13, width = 15)
print(p)
dev.off()
# }}}
# group lasso
# {{{
WD = "~/@patrick/SF9495/integration_analysis/lasso_modeling/glasso/cumulative/stability_threshold"
setwd(WD)
base::load(paste(WD, "workspace_after_bootstrap.Robj", sep = "/"), verbose = TRUE)
base::load(paste(WD, "workspace_after_bootstrapPerm.Robj", sep = "/"), verbose = TRUE)
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
clones = c("Clone.1", "Clone.3", "Clone.4", "Clone.5", "Clone.6")
names(tvalReal) = clones
names(tvalPerm) = clones
tvalRealDf = reshape2::melt(tvalReal)
tvalPermDf = reshape2::melt(tvalPerm)
tvalRealDf = data.table(tvalRealDf, Permuted = rep("Real", nrow(tvalRealDf)))
tvalPermDf = data.table(tvalPermDf, Permuted = rep("Permuted", nrow(tvalPermDf)))
tvalPlot = rbind(tvalRealDf, tvalPermDf)
tvalPlot$L1 = gsub("\\.", " ", tvalPlot$L1)
tvalPlot$Permuted = factor(tvalPlot$Permuted, levels = c("Real", "Permuted"))
tvalPlot$L1 = gsub("Clone.1", "Clone 1", tvalPlot$L1)
tvalPlot$L1 = factor(tvalPlot$L1, levels = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6"))

clones = c("Clone 1", "Clone 3", "Clone 4", "Clone 5", "Clone 6")
pval = lapply(clones, function(x) {
    samp3k = seq(1, length(tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Real"]))
    z = ad.test(tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Real"][samp3k], tvalPlot$value[tvalPlot$L1 == x & tvalPlot$Permuted == "Permuted"][samp3k])$ad[2, 3]
    return(z)
})

p2 = ggplot(tvalPlot, aes(x = L1, y = value, fill = Permuted)) +
    geom_violin() +
    theme_classic() +
    oldham_theme() +
    scale_y_continuous(limits = c(-45, 45), breaks = seq(-45, 45, 15)) +
    labs(x = "", y = "T-values", fill = "", title = "T-values of group lasso model") + # title='T-values of cumulative group-lasso model', subtitle='P-values from two-sided Anderson-Darling test') +
    scale_fill_manual(values = c("black", "white")) +
    theme(
        legend.spacing.y = unit(1, "cm"),
        axis.text.x = element_text(angle = 90)
    ) +
    geom_signif(
        y_position = c(rep(35, 5)),
        xmin = seq(0.8, 4.8, 1),
        xmax = seq(1.2, 5.2, 1),
        annotation = signif(unlist(pval), 1),
        tip_length = 0,
        size = 2,
        textsize = 10,
        vjust = -.6,
        family = "NimbusSan"
    ) +
    guides(fill = guide_legend(byrow = TRUE))

library("gridExtra")
p1 = p1 +
    theme(axis.text.x = element_blank())
pdf("~/@patrick/SF9494/figures/tvalue_lasso.pdf", height = 13 * 2, width = 13)
grid.arrange(arrangeGrob(p1, p2, nrow = 2), ncol = 1)
dev.off()
p1_alt = p1 + theme(legend.position = "none", plot.margin = margin(b = 0)) + scale_y_continuous(limits = c(-45, 45), breaks = seq(-45, 45, 45))
p2_alt = p2 + theme(legend.position = "none", plot.margin = margin(t = 0)) + scale_y_continuous(limits = c(-45, 45), breaks = seq(-45, 45, 45))
pdf("~/@patrick/SF9495/figures/tvalue_lasso_no_leg.pdf", height = 13, width = 13 * .8)
grid.arrange(arrangeGrob(p1_alt, p2_alt, nrow = 2, heights = c(4 / 10, 5 / 10), padding = 0), ncol = 1)
dev.off()
# }}}
