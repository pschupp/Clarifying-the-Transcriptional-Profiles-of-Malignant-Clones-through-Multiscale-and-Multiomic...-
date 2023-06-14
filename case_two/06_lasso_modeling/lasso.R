library('future')
library('future.apply')
library('data.table')
library('reshape2')
library('seagull')
library('ggplot2')
plan(multicore, workers=18)
options(future.globals.maxSize = as.numeric('+Inf'))
setwd('~/@patrick/SF10711/figures/fig6')
options(future.rng.onMisuse='ignore')
# load clonal abundance data
# {{{
barOut=data.frame(fread('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv'))
# construct cumulative clonal abundances by addition
barOutInNonCum=barOut[,-c(1,2)]
colnames(barOutInNonCum)=c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
barOutInNonCum=data.frame(apply(barOutInNonCum, 2, scale))
barOutIn=data.frame(barOut[,-c(1,2,3,4)], Malignant=1-barOut$nonmalignant, Clone.2=(barOut$clone.2+barOut$clone.3))
barOutIn=barOutIn[,c(4,5,1,2,3)]
colnames(barOutIn)=c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
barOutIn=data.frame(apply(barOutIn, 2, scale))
# }}}
# load RNA-seq data
# {{{
rna=fread('/mnt/bdata/@patrick/SF10711/rna.seq/expression_matrices/Normalized_read_counts_using_RUVg_ERCC_K10Factors.csv', drop=seq(2,6))
colnames(rna)[1]='Gene'
meanExpr=apply(rna[,-seq(1,6)], 1, var)
rna=rna[order(meanExpr, decreasing=T),]
rnarSums=apply(rna[,-1], 1, sum)
rna=rna[-which(rnarSums<quantile(rnarSums,0.15)),]
rna=rna[-grep('^ERCC-', rna$Gene),]
rnarSumsN=apply(rna[,-1], 1, sum)
rna=data.frame(Gene=rna$Gene, as.data.frame(rna[,-1]+1))
rna=rna[,c(1, which(colnames(rna) %in% barOut$ampseq.sample))]
rnaScale=t(future_apply(rna[,-c(1)], 1, log2))
# correct line below and make sure that the sections line up as expected
colnames(rnaScale)=gsub('SF10711)', '', colnames(rnaScale))
# }}}
# make correlation heatmaps
# {{{
heatmapOldham = function(df, sublabel){
    corPlot=data.frame(cor(df, method='spearman'))
    corOrder=hclust(as.dist(1-corPlot))$order
    corPlot = corPlot[corOrder, corOrder]
    corOrdN = rownames(corPlot)
    corPlot=melt(data.frame(rownames(corPlot), corPlot))
    colnames(corPlot)=c('Var1', 'Var2', 'value')
    corPlot$Var1 = factor(corPlot$Var1, levels= corOrdN)
    corPlot$Var2 = factor(corPlot$Var2, levels= corOrdN)
    label=signif(corPlot$value, 2)
    print(ggplot(corPlot, aes(x=Var1, y=Var2, fill=value)) +
        theme_classic()+
        oldham_theme()+
        geom_tile() +
        geom_text(aes(Var1, Var2, label=label), color='black',size=7) +
        scale_fill_distiller(palette='RdBu', breaks=seq(-1,1,1), limits=c(-1,1))+
        labs(title='Correlation plot\nof modeling vectors', subtitle=sublabel, x='', y='', fill='Pearson\ncorrelation') +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    )
}
pdf('cor_heatmap.pdf', height=13, width=13)
heatmapOldham(barOutIn, 'Cumulative model')
dev.off()
# }}}
# group/lasso lasso model implemented via seagull
# {{{
# group lasso is probably less applicable in this case as there are fewer genes correlated with it, and there is relatively little variation in it.
library('seagull')
library('metap')
library('progressr')
handlers(global = TRUE)
handlers("progress")
# config options
boots = 100# number of bootstraps
subsetSamp = seq(1,nrow(rnaScale)) # subset which genes to model
coefOut = pvalOut = list()
samp=seq(1,boots)
getAIC = function(i, reglas){ # function that returns a list of the linear model and the akike information criterion (AIC)
    test = reglas[[i]]$random_effects
    testIndeces = unique(apply(test, 1, function(x) which(x != 0)))
    if(any(lapply(testIndeces, length) == 0)){
        testIndeces = testIndeces[-which(lapply(testIndeces, length) == 0)]
    }
    testTemp = lapply(testIndeces, function(x) data.frame(y=(rnaScale[i,]), barOutIn[,x]))
    if(any(unlist(lapply(testIndeces, length),recursive=T)==1)){
        for(ea in which(unlist(lapply(testIndeces, length))==1)){
            colnames(testTemp[[ea]])=c('y', colnames(barOutIn)[unlist(testIndeces[ea])])
        }
    }
    testLm = lapply(testTemp, function(x) lm(y ~ ., data=x))
    testAIC = lapply(testLm, BIC) # get BIC, select the model with the most valuable factors relative to the number of coefficients. AIC most strongly penalizes extra coefficients realtive to BIC or adjusted r-squared
    return(list(testLm, testAIC))
}

bootstrapGlasso = function(samp) {
    p = progressor( along = samp)
    bootOut  = future_lapply(samp, function(samples) {
        sampB=sample(seq(1,ncol(rnaScale)), ncol(rnaScale), replace=T)
        reglas = apply(rnaScale[subsetSamp,sampB], 1, function(y)
                     seagull(y = y,
                             Z = as.matrix(barOutIn[sampB,]),
                             alpha = 1, # exclusive LASSO
                             standardize = T,
                             max_iter = 500)
        )
        out = future_lapply(seq(1, length(reglas)), function(x) getAIC(x, reglas)) # a list for every gene of the regression model and AIC
        aicMin=unlist(future_lapply(seq(1,length(reglas)), function(x) which.min(unlist(out[[x]][[2]]))))
        coefOut = future_lapply(seq(1, length(reglas)), function(i) out[[i]][[1]][[aicMin[i]]]$coefficients[-1])
        tvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$coefficients[-1,3])
        FvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$fstatistic)
        pvalOut = future_lapply(FvalOut, function(x) 1-pf(x[1], x[2], x[3]))
        cat(paste('Finished bootstrap...', samples, '\n', sep=''))
        return(list(coefOut, pvalOut, tvalOut))
    })
    return(bootOut)
}
bootOut = bootstrapGlasso(samp)
outBoot=value(bootOut)
coefOut=lapply(outBoot, function(x) x[[1]])
pvalOut=lapply(outBoot, function(x) x[[2]])
tvalOut=lapply(outBoot, function(x) x[[3]])
# get stability which is the number of times the most frequent combination of factors is chosen
stability=unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) max(table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse='_')))))))
# p-value across all experiments averaged using Fisher's method
pvalBag=unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) sumlog(unlist(lapply(pvalOut, function(x) x[i])))$p))
# average length of chosen model
lengthBag=unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) mean(unlist(lapply(lapply(coefOut, function(x) x[i][[1]]), length)))))
save(list=c('coefOut', 'pvalOut', 'tvalOut', 'stability', 'pvalBag', 'lengthBag', 'rnaScale', 'rna', 'barOutIn'), file='workspace_after_bootstrap.Robj')
# evaluation functions

plotDf=data.frame(as.factor(stability), pvalBag, lengthBag)
pdf('bagged_values.pdf', height=13, width=13)
print(ggplot(plotDf, aes(x=lengthBag))+
        labs(x='Number of terms', y='Count', title='Histogram of number of terms\nin each bootstrap for each gene') +
        geom_histogram() +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDf, aes(x=stability))+
        labs(y='Count', x='Stability', title='Histogram of number of stability\nin each bootstrap for each gene') +
        geom_histogram() +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDf, aes(x=as.factor(stability), y=lengthBag))+
        labs(x='Stability', y='Number of terms', title='Change in average number of\ncoefficients as the stability\nincreases') +
        geom_boxplot(notch=T) +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDf, aes(x=as.factor(stability), y=(log(pvalBag+1E-17)*-1)))+
        labs(x='Stability', y='-log(p-value)', title="P-value as stability increases")+
        geom_boxplot(notch=T) +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDf, aes(x=as.factor(round(lengthBag)), y=(log(pvalBag+1E-17)*-1)))+
        labs(x='Rounded number of terms', y='-log(p-value)', title='P-value as the number\nof terms increases') +
        geom_boxplot(notch=T) +
        theme_classic() +
        oldham_theme())
dev.off()
# now need to do permuted version for control
# will permute the clonal vectors barOutIn each time prior to running the model 
bootstrapGlasso = function(samp) {
    p = progressor( along = samp)
    bootOut= future_lapply(samp, function(samples) {
        sampA=sample(seq(1,ncol(rnaScale)), ncol(rnaScale), replace=T)
        sampB=sample(seq(1,ncol(rnaScale)), ncol(rnaScale), replace=T)
        reglas = apply(rnaScale[subsetSamp,sampB], 1, function(y)
                     seagull(y = y,
                             Z = as.matrix(barOutIn[sampA,]),
                             alpha = 1, # exclusive LASSO
                             standardize = T,
                             max_iter = 500)
        )
        out= future_lapply(seq(1, length(reglas)), function(x) getAIC(x, reglas)) # a list for every gene of the regression model and AIC
        aicMin=unlist(future_lapply(seq(1,length(reglas)), function(x) which.min(unlist(out[[x]][[2]]))))
        coefOut = future_lapply(seq(1, length(reglas)), function(i) out[[i]][[1]][[aicMin[i]]]$coefficients[-1])
        tvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$coefficients[-1,3])
        FvalOut = future_lapply(seq(1, length(reglas)), function(i) summary(out[[i]][[1]][[aicMin[i]]])$fstatistic)
        pvalOut = future_lapply(FvalOut, function(x) 1-pf(x[1], x[2], x[3]))
        cat(paste('Finished bootstrap...', samples, '\n', sep=''))
        return(list(coefOut, pvalOut, tvalOut))
    })
    return(bootOut)
}
bootOutPerm = bootstrapGlasso(samp)
outPermBoot=value(bootOutPerm)
coefOutPerm=lapply(outPermBoot, function(x) x[[1]])
pvalOutPerm=lapply(outPermBoot, function(x) x[[2]])
tvalOutPerm=lapply(outPermBoot, function(x) x[[3]])
# get stability which is the number of times the most frequent combination of factors is chosen
stabilityPerm=unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) max(table(unlist(lapply(coefOutPerm, function(x) paste(names(x[i][[1]]), collapse='_')))))))
# p-value across all bootstraps for each gene averaged using Fisher's method
pvalBagPerm=unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) sumlog(unlist(lapply(pvalOutPerm, function(x) x[i])))$p))
# average length of chosen model for each gene
lengthBagPerm=unlist(future_lapply(seq(1, nrow(rnaScale)), function(i) mean(unlist(lapply(lapply(coefOutPerm, function(x) x[i][[1]]), length)))))
save(list=c('coefOutPerm', 'pvalOutPerm', 'tvalOutPerm', 'stabilityPerm', 'pvalBagPerm', 'lengthBagPerm', 'rnaScale', 'rna', 'barOutIn'), file='workspace_after_bootstrapPerm.Robj')
plotDfPerm=data.frame(as.factor(stabilityPerm), pvalBagPerm, lengthBagPerm)
pdf('bagged_values_perm.pdf', height=13, width=13)
print(ggplot(plotDfPerm, aes(x=lengthBagPerm))+
        labs(x='Number of terms', y='Count', title='Histogram of number of terms\nin each bootstrap for each gene') +
        geom_histogram() +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDfPerm, aes(x=stabilityPerm))+
        labs(y='Count', x='Stability', title='Histogram of number of stability\nin each bootstrap for each gene') +
        geom_histogram() +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDfPerm, aes(x=as.factor(stabilityPerm), y=lengthBagPerm))+
        labs(x='Stability', y='Number of terms', title='Change in average number of\ncoefficients as the stability\nincreases') +
        geom_boxplot(notch=T) +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDfPerm, aes(x=as.factor(stabilityPerm), y=(log(pvalBagPerm+1E-17)*-1)))+
        labs(x='Stability', y='-log(p-value)', title="P-value as stability increases")+
        geom_boxplot(notch=T) +
        theme_classic() +
        oldham_theme())
print(ggplot(plotDfPerm, aes(x=as.factor(round(lengthBagPerm)), y=(log(pvalBagPerm+1E-17)*-1)))+
        labs(x='Rounded number of terms', y='-log(p-value)', title='P-value as the number\nof terms increases') +
        geom_boxplot(notch=T) +
        theme_classic() +
        oldham_theme())
dev.off()
# }}}
