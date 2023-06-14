library("future")
library("future.apply")
library("data.table")
library("reshape2")
library("seagull")
library("metap")
library("progressr")
library("progress")
library("ggplot2")
plan(multicore, workers = 15)
options(future.globals.maxSize = as.numeric("+Inf"))
setwd("/mnt/bdata/@patrick/SF9495/integration_analysis/glasso")
options(future.rng.onMisuse = "ignore")
handlers(global = TRUE)
handlers("progress")

# import RNA-seq data
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
rnaScale = t(future_apply(rna[, -c(seq(1, 6), 69, 70)], 1, scale))
colnames(rnaScale) = gsub("SF9495_", "", colnames(rnaScale))
# }}}

# import clonal abundance
# {{{{
barOut = fread("~/@patrick/SF9495/figures/sf9495_cellular_abundance.csv")
barOut = data.frame(barOut, Malignant = 1 - barOut$Nonmalignant, clone2Add = (barOut$clone.2 + barOut$clone.3 + barOut$clone.4 + barOut$clone.5 + barOut$clone.6), clone3Add = (barOut$clone.3 + barOut$clone.4 + barOut$clone.5))
colnames(barOut)[1:9] = c("section", "Nonmalignant", "clone 1", "clone 2", "clone 3", "clone 4", "clone 5", "clone 6", "Malignant")
barOut = data.frame(barOut)
barOut[, seq(2, ncol(barOut))] = apply(barOut[, seq(2, ncol(barOut))], 2, scale)
barOutIn = data.frame(barOut$Malignant, barOut$clone2Add, barOut$clone3Add, barOut$clone.5, barOut$clone.6, barOut$clone.4)
colnames(barOutIn) = gsub("barOut\\.", "", colnames(barOutIn))
colnames(barOutIn) = c("Malignant", "Clone.2", "Clone.3", "Clone.5", "Clone.6", "Clone.4")
evoOrder = c(1, 2, 3, 6, 4, 5)

cnvDf = data.frame(chr = c(10, 2, 2), start = c(93816, 232035429, 46043), stop = c(39047679, 243072740, 10956969), type = c("amplification", "deletion", "deletion"), clone = c("clone 5", "clone 6", "clone 3,4,5"), clone_ind = c(4, 5, "3,4,6"), name = c("Chr10", "Chr2q", "chr2p"))
genes = as.character(unlist(rna[, 2]))
# }}}

# config options
boots = 100 # number of bootstraps
subsetSamp = seq(1, nrow(rnaScale)) # subset which genes to model
# perform bootstrap lasso
# {{{
coefOut = pvalOut = list()
samp = seq(1, boots)
getAIC = function(i, reglas) { # function that returns a list of the linear model and the akike information criterion (AIC)
    test = reglas[[i]]$random_effects
    testIndeces = unique(apply(test, 1, function(x) which(x != 0)))
    if (any(lapply(testIndeces, length) == 0)) {
        testIndeces = testIndeces[-which(lapply(testIndeces, length) == 0)]
    }
    testTemp = lapply(testIndeces, function(x) data.frame(y = (rnaScale[i, ]), barOutIn[, x]))
    if (any(unlist(lapply(testIndeces, length), recursive = T) == 1)) {
        for (ea in which(unlist(lapply(testIndeces, length)) == 1)) {
            colnames(testTemp[[ea]]) = c("y", colnames(barOutIn)[unlist(testIndeces[ea])])
        }
    }
    testLm = lapply(testTemp, function(x) lm(y ~ ., data = x))
    testAIC = lapply(testLm, BIC) # get BIC, select the model with the most valuable factors relative to the number of coefficients. AIC most strongly penalizes extra coefficients realtive to BIC or adjusted r-squared
    return(list(testLm, testAIC))
}

bootstrapGlasso = function(samp) {
    p = progressor(along = samp)
    bootOut = future_lapply(samp, function(samples) {
        sampB = sample(seq(1, ncol(rnaScale)), ncol(rnaScale), replace = T)
        reglas = apply(rnaScale[subsetSamp, sampB], 1, function(y) {
            seagull(
                y = y,
                Z = as.matrix(barOutIn[sampB, ]),
                groups = c(1, rep(2, 5)),
                alpha = 1, # exclusive LASSO
                standardize = F,
                max_iter = 500
            )
        })
        out = future_lapply(seq(1, length(reglas)), function(x) getAIC(x, reglas)) # a list for every gene of the regression model and AIC
        aicMin = unlist(future_lapply(seq(1, length(reglas)), function(x) which.min(unlist(out[[x]][[2]]))))
        coefOut = future_lapply(seq(1, length(reglas)), function(i) out[[i]][[1]][[aicMin[i]]]$coefficients[-1])
        tvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$coefficients[-1, 3])
        FvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$fstatistic)
        pvalOut = future_lapply(FvalOut, function(x) 1 - pf(x[1], x[2], x[3]))
        cat(paste("Finished bootstrap...", samples, "\n", sep = ""))
        return(list(coefOut, pvalOut, tvalOut))
    })
    return(bootOut)
}
bootOut = bootstrapGlasso(samp)
outBoot = value(bootOut)
coefOut = lapply(outBoot, function(x) x[[1]])
pvalOut = lapply(outBoot, function(x) x[[2]])
tvalOut = lapply(outBoot, function(x) x[[3]])
# get stability which is the number of times the most frequent combination of factors is chosen
stability = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) max(table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse = "_")))))))
coefNames = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse = "_"))))))
tN = sort(table(names(coefNames)), decreasing = T)
highStab = stability > 81
coefNamesH = unlist(future_lapply(which(highStab), function(i) table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse = "_"))))))
tH = sort(table(names(coefNamesH)), decreasing = T)
lowStab = stability < 36
coefNamesL = unlist(future_lapply(which(lowStab), function(i) table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse = "_"))))))
tL = sort(table(names(coefNamesL)), decreasing = T)
HhighStab = stability > 72
coefNamesHH = unlist(future_lapply(which(HhighStab), function(i) table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse = "_"))))))
tHH = sort(table(names(coefNamesH)), decreasing = T)

# p-value across all experiments averaged using Fisher's method
pvalBag = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) sumlog(unlist(lapply(pvalOut, function(x) x[i])))$p))
# average length of chosen model
lengthBag = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) mean(unlist(lapply(lapply(coefOut, function(x) x[i][[1]]), length)))))
save(list = c("coefOut", "pvalOut", "tvalOut", "stability", "pvalBag", "lengthBag", "rnaScale", "rna", "barOutIn"), file = "workspace_after_bootstrap.Robj")
# evaluation functions

plotDf = data.frame(as.factor(stability), pvalBag, lengthBag)
pdf("bagged_values.pdf", height = 13, width = 13)
print(ggplot(plotDf, aes(x = lengthBag)) +
    labs(x = "Number of terms", y = "Count", title = "Histogram of number of terms\nin each bootstrap for each gene") +
    geom_histogram() +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDf, aes(x = stability)) +
    labs(y = "Count", x = "Stability", title = "Histogram of number of stability\nin each bootstrap for each gene") +
    geom_histogram() +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDf, aes(x = as.factor(stability), y = lengthBag)) +
    labs(x = "Stability", y = "Number of terms", title = "Change in average number of\ncoefficients as the stability\nincreases") +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDf, aes(x = as.factor(stability), y = (log(pvalBag + 1E-17) * -1))) +
    labs(x = "Stability", y = "-log(p-value)", title = "P-value as stability increases") +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDf, aes(x = as.factor(round(lengthBag)), y = (log(pvalBag + 1E-17) * -1))) +
    labs(x = "Rounded number of terms", y = "-log(p-value)", title = "P-value as the number\nof terms increases") +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme())
dev.off()

# now need to do permuted version for control
# will permute the clonal vectors barOutIn each time prior to running the model
bootstrapGlasso = function(samp) {
    p = progressor(along = samp)
    bootOut = future_lapply(samp, function(samples) {
        sampA = sample(seq(1, ncol(rnaScale)), ncol(rnaScale), replace = T)
        sampB = sample(seq(1, ncol(rnaScale)), ncol(rnaScale), replace = T)
        reglas = apply(rnaScale[subsetSamp, sampB], 1, function(y) {
            seagull(
                y = y,
                Z = as.matrix(barOutIn[sampB, ]),
                groups = c(1, rep(2, 5)),
                alpha = 1, # exclusive LASSO
                standardize = F,
                max_iter = 500
            )
        })
        out = future_lapply(seq(1, length(reglas)), function(x) getAIC(x, reglas)) # a list for every gene of the regression model and AIC
        aicMin = unlist(future_lapply(seq(1, length(reglas)), function(x) which.min(unlist(out[[x]][[2]]))))
        coefOut = future_lapply(seq(1, length(reglas)), function(i) out[[i]][[1]][[aicMin[i]]]$coefficients[-1])
        tvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$coefficients[-1, 3])
        FvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$fstatistic)
        pvalOut = future_lapply(FvalOut, function(x) 1 - pf(x[1], x[2], x[3]))
        cat(paste("Finished bootstrap...", samples, "\n", sep = ""))
        return(list(coefOut, pvalOut, tvalOut))
    })
    return(bootOut)
}
bootOutPerm = bootstrapGlasso(samp)
outPermBoot = value(bootOutPerm)
coefOutPerm = lapply(outPermBoot, function(x) x[[1]])
pvalOutPerm = lapply(outPermBoot, function(x) x[[2]])
tvalOutPerm = lapply(outPermBoot, function(x) x[[3]])
# get stability which is the number of times the most frequent combination of factors is chosen
stabilityPerm = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) max(table(unlist(lapply(coefOutPerm, function(x) paste(names(x[i][[1]]), collapse = "_")))))))
# p-value across all bootstraps for each gene averaged using Fisher's method
pvalBagPerm = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) sumlog(unlist(lapply(pvalOutPerm, function(x) x[i])))$p))
# average length of chosen model for each gene
lengthBagPerm = unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) mean(unlist(lapply(lapply(coefOutPerm, function(x) x[i][[1]]), length)))))
save(list = c("coefOutPerm", "pvalOutPerm", "tvalOutPerm", "stabilityPerm", "pvalBagPerm", "lengthBagPerm", "rnaScale", "rna", "barOutIn"), file = "workspace_after_bootstrapPerm.Robj")

plotDfPerm = data.frame(as.factor(stabilityPerm), pvalBagPerm, lengthBagPerm)
pdf("bagged_values_perm.pdf", height = 13, width = 13)
print(ggplot(plotDfPerm, aes(x = lengthBagPerm)) +
    labs(x = "Number of terms", y = "Count", title = "Histogram of number of terms\nin each bootstrap for each gene") +
    geom_histogram() +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDfPerm, aes(x = stabilityPerm)) +
    labs(y = "Count", x = "Stability", title = "Histogram of number of stability\nin each bootstrap for each gene") +
    geom_histogram() +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDfPerm, aes(x = as.factor(stabilityPerm), y = lengthBagPerm)) +
    labs(x = "Stability", y = "Number of terms", title = "Change in average number of\ncoefficients as the stability\nincreases") +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDfPerm, aes(x = as.factor(stabilityPerm), y = (log(pvalBagPerm + 1E-17) * -1))) +
    labs(x = "Stability", y = "-log(p-value)", title = "P-value as stability increases") +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme())
print(ggplot(plotDfPerm, aes(x = as.factor(round(lengthBagPerm)), y = (log(pvalBagPerm + 1E-17) * -1))) +
    labs(x = "Rounded number of terms", y = "-log(p-value)", title = "P-value as the number\nof terms increases") +
    geom_boxplot(notch = T) +
    theme_classic() +
    oldham_theme())
dev.off()
# }}}
