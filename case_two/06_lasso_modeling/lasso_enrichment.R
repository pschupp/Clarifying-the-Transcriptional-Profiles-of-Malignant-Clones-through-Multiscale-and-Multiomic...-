# load in R workspace of lasso
# {{{
base::load('~/@patrick/SF10711/figures/fig6/workspace_after_bootstrap.Robj')
base::load('~/@patrick/SF10711/figures/fig6/workspace_after_bootstrapPerm.Robj')
setwd('~/@patrick/SF10711/figures/fig6/')
# }}}
# coefOut - length 100, for all bootstraps coefficients for model with minimum AIC, this model is the relaxed lasso model, so no restriction on scale of coefficients
# pvalOut - length 100, for all bootstraps model p-value for model with minimum AIC
# tvalsOut - for all bootstraps, t-values for all coefficients in the model, length will the same for each member as coefOut
# lengthBag - average number of coefficients for each gene across all bootstraps
# stability - highest indcidence of a combination of clone vector (independent variable) across all bootsraps
# pvalBag - p-value across all bootstraps for each gene averaged using Fisher's method
# dependencies
# {{{
library('future')
library('future.apply')
library('ggplot2')
library('qvalue')
plan(multicore, workers=18)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
# mode (most prevalent) factor combination for the model
bootMode=future_lapply(seq(1, nrow(rnaScale)), function(i) names(which.max(table(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse='_')))))))
# }}}
# stability threshold is greater than 44
# get genes with +/- correlations
# {{{
# first, for each gene, select only the bootstraps which are the boot mode
modeSelec=future_lapply(seq(1, nrow(rnaScale)), function(i) which(unlist(lapply(coefOut, function(x) paste(names(x[i][[1]]), collapse='_')))==bootMode[i]))
# now need to use modeSelec to only sample from the correct bootstraps
future_lapply(seq(1,nrow(rnaScale)), function(i) x[[1]])
pvalBagMode=unlist(future_lapply(seq(1,nrow(rnaScale)), function(i) mean(unlist(lapply(pvalOut[modeSelec[[i]]], function(x) x[[i]])))))
pvalBagNoMode=unlist(future_lapply(seq(1,nrow(rnaScale)), function(i) mean(unlist(lapply(pvalOut[-modeSelec[[i]]], function(x) x[[i]])))))
stabThresh=44
x=sort(table(unlist(bootMode[stability>stabThresh & qvalue(pvalBag)$qvalue<0.05])),decreasing=T)
# get list of genes with positive and negative coefficients
bootGenes=list()
i=1
j=1
cloneVec=sort(unique(unlist(strsplit(unique(unlist(bootMode)),'_'))))
for(cloneName in cloneVec){
    temp1=which(bootMode[stability>stabThresh & qvalue(pvalBag)$qvalue<0.05]==cloneName)
    temp2=seq(1,nrow(rnaScale))[stability>stabThresh & qvalue(pvalBag)$qvalue<0.05][temp1]
    temp3=apply(rnaScale[temp2,], 1, function(rnaRow) lm(rnaRow ~ barOutIn[,which(colnames(barOutIn) == cloneName)]))
    temp4=unlist(lapply(temp3, function(x) x$coefficients[2]))
    t4pos=which(temp4>0)
    t4neg=which(temp4<0)
    genesP=rna$Gene[temp2][t4pos]
    genesN=rna$Gene[temp2][t4neg]
    bootGenes[[j]]=genesP
    bootGenes[[j+1]]=genesN
    j=j+2
    i=i+1
}
# }}}
# fisher's exact test enrichment calculation of oldham and broad sets
# {{{
names(bootGenes)=paste(rep(cloneVec,each=2), rep(c('_positive', '_negative'), times=length(cloneVec)), sep='')
source('~/code/git/GSEA_generic/enrichment_functions.R')
x=enrichment_man(enrich_vec=bootGenes, all_genes=rna$Gene, enrichDir='/home/shared/genesets/genesets_slim')
write.csv(x, file='enrich_results_boot_model.csv', row.names=F)
library('GSEABase')
broadSets=getBroadSets('/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml')
y=enrichment_man_broad(bootGenes, rna$Gene, broadSets=broadSets)
write.csv(y, file='enrich_results_broad_boot_model.csv', row.names=F)
# }}}
# Fisher's exact test enrichment result of CNV sets
# {{{
library('data.table')
gtf=fread('/home/shared/hg_align_db/GRCh37.38/gencode.v38lift37.annotation.gene_only.gtf')
gtf=gtf[,-c(2,3,6,8)]
colnames(gtf)=c('chr', 'start', 'end','strand', 'gene')
gtf$gene=gsub('";.*', '', gsub('.*gene_name "', '', gtf$gene))

cnvDf=data.frame(chr=c(2,18, 11, 7, 8, 9),
                    start=c(98204464, 41955206, 55838392, 192969, 11977, 12705),
                    stop=c(224844864, 77927112, 134605703, 159144957, 146289882, 141124247),
                    type=c('deletion', 'deletion', 'amplification', 'amplification', 'amplification', 'amplification'),
                    clone=c('Clone 3', 'Clone 3', 'Clone 5', 'Clone 1', 'Clone 1', 'Clone 1'),
                    clone_ind=c(3, 3, 5, 1, 1, 1),
                    names=c('Chr2q', 'Chr18', 'Chr11q', 'Chr7', 'Chr8', 'Chr9'))
cnvGenes=list()
for(i in seq(1,nrow(cnvDf))){
    cnvGenes[[i]]=gtf$gene[which(gtf$chr==paste0('chr', cnvDf$chr[i]) & gtf$start>cnvDf$start[i] & gtf$end<cnvDf$stop[i])]
}
names(cnvGenes)=c('Chr2q deletion', 'Chr18 deletion', 'Chr11q amplification', 'Chr7 amplification', 'Chr8 amplification', 'Chr9 amplification')
cnvSets=data.frame(gs_name=rep(names(cnvGenes), unlist(lapply(cnvGenes, length))), entrez_gene=unlist(cnvGenes))
# get the cnvGS genes into files so we can run enrichment_man
setwd('~/@patrick/SF10711/figures/fig6/glasso/cnv_gene_set')
for(i in seq_along(cnvGenes)){
    write.table(cnvGenes[[i]], paste0(names(cnvGenes)[i], '.csv'), row.names=F, quote=F, col.names=F, sep=',')
}
names(bootGenes)=paste(rep(cloneVec,each=2), rep(c('_positive', '_negative'), times=length(cloneVec)), sep='')
source('~/code/git/GSEA_generic/enrichment_functions.R')
cnvEnrich=enrichment_man(bootGenes, rna$Gene, '~/@patrick/SF10711/figures/fig6/glasso/cnv_gene_set')
setwd('~/@patrick/SF10711/figures/fig6/glasso')
write.table(cnvEnrich, file='cnvEnrichment_results_boot.csv', row.names=F, sep=',', quote=F)
# }}}
