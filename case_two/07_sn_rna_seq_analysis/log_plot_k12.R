# environment setup
# {{{
library('data.table')
library('ggplot2')
library('future')
library('future.apply')
library('qvalue')
library('dendextend')
library('RColorBrewer')
library('grid')
library('gridExtra')
library('limma')
library('ComplexHeatmap')
library('parallelDist')
library('kSamples')
library('ggsignif')
library('circlize')
library('phylogram')
library('slingshot')
library('igraph')
plan(multicore, workers=20)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
WD = '~/@patrick/SF10711/sn.rna.seq/figures/'
# }}}
# panel A - violin plots showing re distribution of t-values for the geneset versus versus all other clones 
# {{{
# get FM clone genesets and nonmalignant genesets
# {{{
clones = c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
fmGenes = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
clones=c('red', 'violet', 'black', 'ivory', 'lightcyan')
fmGenes = lapply(clones, function(clone) fmGenes$Gene[which(fmGenes$'TopModPosFDR_0.12' == clone)])
names(fmGenes) = clones

fidGenes = lapply(list.files('/home/shared/genesets/fidelity_cell_types/hifi200', full.names=T), function(x) as.character(fread(x)[[1]]))
names(fidGenes) = gsub('_.*', '', list.files('/home/shared/genesets/fidelity_cell_types/hifi200'))
fidGenes= fidGenes[c(1, 10, 7, 9, 5)]
# we can instead get the genes from the dataset
clones = c('Astrocyte', 'Oligodendrocyte', 'Microglia', 'Neuron', 'Endothelial cell')
bulkGenes = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
clones=c('purple', 'tan', 'yellow', 'green', 'orange')
bulkGenes = lapply(clones, function(clone) bulkGenes$Gene[which(bulkGenes$'TopModPosBC_3.63e-08' == clone)])
names(bulkGenes) = clones
set1 = c(fmGenes, bulkGenes)
# }}}
# get expression matrix and cluster assignments
# {{{
dist = as.numeric(unlist(fread('~/@patrick/SF10711/sn.rna.seq/sanity/cell_cell_distance_with_errorbar_s2n_gt_1.txt')))
# read in expression matrix to determine wanted dimensions
expr = data.frame(fread('~/@patrick/SF10711/sn.rna.seq/sanity/log_transcription_quotients.txt'))
expr1 =expr[,1]
expr = expr[,-1]
expr=data.table(expr1, apply(expr, 2, function(x) exp(x))*1E6)
# format into distance matrix
zeros = seq(1,ncol(expr)-1)
delta = seq((ncol(expr)-3), 1)
indeces = c(1)
for(i in seq(1, ncol(expr)-3)){
    if(i==1){
        indeces [i+1] = indeces[i] + delta[i]
    } else{
        indeces [i+1] = indeces[i] + delta[i] +1
    }
}
indeces[1]=0
indeces[length(indeces)+1]=indeces[length(indeces)]+1
distMat=lapply(seq_along(indeces[-length(indeces)]), function(i) {
    out=dist[(indeces[i]+1) : indeces[i+1]]
    if(zeros[i] != 809){
        out=unlist(c(rep(0, zeros[i]),out))
    } else{
        out=unlist(c(rep(0, zeros[i]-1), out[2]))
    }
    return(out)
})
distMat=data.frame(do.call(cbind, distMat))
distMat=cbind(distMat, rep(0, nrow(distMat)))
# now need to mirror distance matrix
for(i in seq_len(nrow(distMat))) {
    for(j in seq_len(ncol(distMat))) {
        distMat[i, j] = distMat[j, i]
    }
}
distDf = as.dist(distMat)
# }}}
# harmonize gene names
# {{{
gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf', skip=5)
colnames(gtf)=c("chr", "db", "type", "start", "end", "score", "strand", "frame", "extra")
gtf=gtf[which(gtf$type=="gene"),]
gtf.ENSG=gsub(".*gene_id\\s\"\\\"*|\\..*", "", gtf$extra)
gtf.type=gsub(".*gene_type\\s\"\\\"*|\".*", "", gtf$extra)
gtf.name=gsub(".*gene_name\\s\"\\\"*|\".*", "", gtf$extra)
gtfOut=data.frame(gtf[,1], gtf.ENSG, gtf.type, gtf.name)
gtfOut$gtf.name= alias2SymbolTable(gtfOut$gtf.name, species = 'Hs')
gtfOut$gtf.name[is.na(gtfOut$gtf.name)] = gtfOut$gtf.ENSG[is.na(gtfOut$gtf.name)]
fmGenes=lapply(fmGenes, alias2SymbolTable, species='Hs')

expr1Gene = expr$expr1
expr = data.frame(gtfOut$gtf.name[match(as.character(unlist(expr$expr1)), gtfOut$gtf.ENSG)], expr)
colnames(expr)[1:2] = c('Gene', 'ENSG')
expr$Gene[is.na(expr$Gene)] = expr$ENSG[is.na(expr$Gene)]
exprC = fread('~/@patrick/SF10711/sn.rna.seq/190809_SF013_oldham_july_1/08_expression_matrix/expr.ensg.counts.csv')
exprC = exprC[match(expr$ENSG, exprC$Gene),]
exprC = data.frame(expr$Gene, exprC)
colnames(exprC)[1:2] = c('Gene', 'ENSG')

countThresh = future_apply(exprC[,-c(1,2)], 1, function(x) length(which(x<0.5)) / (ncol(exprC)-2))
threshInd = which(countThresh<0.9)
exprCD = exprC[threshInd,-c(1,2)]
exprD = expr[threshInd,-c(1,2)]
exprSum = future_apply(exprCD, 1, sum)
exprGene = expr$Gene[threshInd]
fmGenes = lapply(fmGenes, function(x) return(x[!(is.na(x))]))
# }}}
# get cluster definitions from hierarchical clustering
# {{{
attr(distDf, 'Labels') =paste0('x', seq(1, length(attr(distDf, 'Labels'))))
distHclust = hclust(distDf, method = 'ward.D')
dendP = as.dendrogram(distHclust)
clust = cutree(distHclust, k=12)
names(clust) = colnames(expr)[-c(1,2)]
# }}}
# Then reuse code to get distributions of t-values
# {{{
tValues = function(exprMat, inds){
    apply(exprMat, 1, function(i) t.test(as.numeric(i[inds]), as.numeric(i[!inds]), paired=F)$statistic)
}
names(set1) = c('Red - clone 1', 'Violet - clone 2', 'Black - clone 3', 'Ivory - clone 4', 'Lightcyan - clone 5', 'Purple - astrocytes', 'Tan - oligodendrocytes', 'Yellow - microglia', 'Green - neurons', 'Orange - endothelial cells')
tVal = lapply(set1, function(set){
    ind2 = exprGene %in% set
    future_lapply(unique(clust), function(clu){
        ind = (clust == clu)
        print(clu)
        list(tValues(exprD[ind2,], inds = ind), tValues(exprD[!ind2,], inds = ind))
    })
})
sampN = 1000
adOut = future_lapply(tVal, function(y){
    future_lapply(y, function(x){
        wilcox.test(x[[1]], x[[2]][sample(seq(1, length(x[[2]])), 100)], exact = F, alternative = 'greater', correct = F)$p.value*10*8
    })
})
adOut = do.call(rbind , adOut)
colnames(adOut) = paste('Cluster', seq_len(ncol(adOut)))
# }}}
# plot
# {{{
ordVec = c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10) 
sigVec = list(2, 1, 7, c(3,5), 10, 9, c(11,8), 4, 12, 6)[ordVec]
annoVec = lapply(sigVec, function(x) rep('***', length(x)))
colTitle = c('#ff7f00','#e6ab02','#999999','#e7298a','#dd95e8', '#1b9e77','#4daf4a','#701d02','#a65729','#1d4a4d')

plotZ = lapply(seq_along(ordVec), function(h){
    i = ordVec[h]
    try = lapply(seq_along(tVal[[1]]), function(j){
        plot = tVal[[i]][[j]]
        plot2 = reshape2::melt(plot)
        plot2$L1 = gsub(1, 'Gene set', plot2$L1)
        plot2$L1 = gsub(2, 'All other genes', plot2$L1)
        plot2$Clone = j
        colnames(plot2) = c('T_value', 'Set', 'Clone')
        plot2
    })
    try = do.call(rbind, try)
    try$Set = factor(try$Set, levels = c('Gene set', 'All other genes'))
    try$Clone= factor(try$Clone, levels = seq(1,length(tVal[[1]])))
    try$T_value[try$T_value>10] = NA
    try$T_value[try$T_value<(-10)] = NA
    ggplot(try, aes(x = Clone, y = T_value, fill = Set))+
        geom_violin(scale='area') +
        scale_y_continuous(limits=c(-15, 17), breaks = seq(-10, 10, 10))+
        labs(y = '', x = '', title = names(tVal)[i], fill='') +
        geom_hline(yintercept=0, color='grey50', size=0.6)+
        scale_fill_manual(values=c('white', 'black')) +
        geom_signif(annotations = annoVec[[h]],
            xmin = sigVec[[h]]-.3, xmax = sigVec[[h]]+.3,  y_position = 11, tip_length = 0, textsize = 15, vjust=.5, size = 1, family = 'NimbusSan'
        ) +
        theme_classic() +
        oldham_theme() +
        theme(axis.text.x = element_text(size=40, color='black',  family='NimbusSan', margin=margin(t=10)),
                axis.text.y = element_text(size=40, color='black',  family='NimbusSan', margin=margin(r = 15)),
                panel.grid.major.x = element_line(colour = "grey50", size = 0.6),
                plot.title = element_text(colour = colTitle[h]),
                plot.margin = unit(c(0,0,0,0), "cm"), 
                legend.position = 'none')
})
pdf(paste0(WD, 'violin_enrich_no_legend.pdf'), height =, 13*1.4, width = 13*2.8)
plotL=list(grobs=plotZ, ncol=2, 
            bottom = textGrob('snRNA-seq clusters',  gp = gpar(fontsize = 65, fontface='bold', fontfamily='NimbusSan')),
            left = textGrob('T-values',  gp = gpar(fontsize = 65, fontface='bold', fontfamily='NimbusSan'), rot = 90)
)
margin = theme()
grid.arrange(do.call(arrangeGrob, plotL), padding = unit(c(0,0,0,0), "cm"))
dev.off()
# }}}
# }}}
# panel B - two scatter plots on top of each other showing that there is alignment of clonal abuandances from bulk and single-nucleus approaches of malignant or non-malignant populations (the latter will be scaled ME for the relevant cell type)
# {{{
# malignant
# {{{
secs = list(c('N701','N702', 'N709'), c('N712', 'N703', 'N704'), c('N705', 'N706', 'N710'), c('N707', 'N708', 'N711'))
secSum = lapply(secs, function(sec) length(grep(paste(sec, collapse = '|'), names(clust))))
clonesInd = c(2, 1, 7, 5, 3, 10)
outP = lapply(seq_along(clonesInd), function(j){
    i = clonesInd[j]
    temp = clust[clust==i]
    temp2 = lapply(seq_along(secs), function(k){
        sec=secs[[k]]
        (length(which(gsub('.*\\.', '', names(temp)) %in% sec))/as.numeric(secSum[k]))*100
    })
    data.frame(clone = rep(i, length(temp2)), percentage = unlist(temp2), section = c('17', '53', '93', '117'))
})
outP = do.call(rbind, outP)
outP = data.frame(source =rep('single_cell', nrow(outP)), outP)
outP$clone = gsub('^2$', 'clone.1', outP$clone)
outP$clone = gsub('^1$', 'clone.2', outP$clone)
outP$clone = gsub('^7$', 'clone.3', outP$clone)
outP$clone = gsub('^5$', 'clone.4', outP$clone)
outP$clone = gsub('^10$', 'clone.5', outP$clone)
freq = read.csv('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv')
freq$ampseq.sample = as.numeric(gsub('.*_' , '', freq$ampseq.sample))
snSlices = c(17, 53, 93, 117)
inds = lapply(snSlices, function(slice) {
    order(abs(slice - freq$ampseq.sample), decreasing=F)[1:4]
})
weights = c(
    length(grep("N701|N702|N709", colnames(expr))),
    length(grep("N712|N703|N704", colnames(expr))),
    length(grep("N705|N706|N710", colnames(expr))),
    length(grep("N707|N708|N711", colnames(expr)))
)
out = do.call(rbind, lapply(seq_along(inds), function(i) {
    x=inds[[i]]
    apply(freq[x, ], 2, mean)
}))
out = rbind(out, apply(out, 2, function(x) {
    mean(
        c(rep(x[1], weights[1]),
        rep(x[2], weights[2]),
        rep(x[3], weights[3]),
        rep(x[4], weights[4]))
        )
    })
)
out = reshape2::melt(out[-5,])
out$Var1 = gsub('^1$','17', out$Var1)
out$Var1 = gsub('^2$','53', out$Var1)
out$Var1 = gsub('^3$','93', out$Var1)
out$Var1 = gsub('^4$','117', out$Var1)
out = out[-grep('sample|nonmalignant', out$Var2),]
out = data.frame(source=rep('bulk', nrow(out)),clone = out$Var2, percentage=out$value*100, section=out$Var1)
out$percentage[3] = out$percentage[3] + 1
outF = rbind(outP, out)
outF$section = factor(outF$section, levels=c('17', '53', '93', '117'))
outF$clone = gsub('clone\\.', 'Clone ', outF$clone)
outF$source= gsub('bulk', 'Bulk', outF$source)
outF$source= gsub('single_cell', 'Single-nucleus', outF$source)
rownames(outF) = seq(1,nrow(outF))
cor.test(outF$percentage[1:4], outF$percentage[21:24])
cor.test(outF$percentage[5:8], outF$percentage[25:28])
cor.test(outF$percentage[9:12], outF$percentage[29:32])
cor.test(outF$percentage[13:16], outF$percentage[33:36])
cor.test(outF$percentage[17:20], outF$percentage[37:40])
outF$clone = gsub('Clone 1', 'Clone 1, R=0.97', outF$clone)
outF$clone = gsub('Clone 2', 'Clone 2, R=0.95', outF$clone)
outF$clone = gsub('Clone 3', 'Clone 3, R=0.94', outF$clone)
outF$clone = gsub('Clone 4', 'Clone 4:1+2, R=0.96', outF$clone)
outF$clone = gsub('Clone 5', 'Clone 5, R=0.94', outF$clone)
pdf(paste0(WD, 'frequency_scatter.pdf'), height=13, width=13*1.2)
ggplot(outF, aes(x =section, y=percentage, group = interaction(clone, source), color=clone, shape=source))+
    geom_point(size=8) +
    geom_line(size=3) +
    labs(x='Section', y='Percentage', title='Bulk and single-nucleus clonal\nfrequency estimates coincide', shape = 'Source of estimate', color = 'Clone') +
    scale_y_continuous(limit=c(0,60), breaks = seq(0, 60, 20)) +
    scale_x_discrete(expand = c(0.01,0.01)) +
    scale_color_manual(values = c('#ff7f00','#999999', '#dd95e8', '#4daf4a', '#a65729'))+
    theme_classic() +
    oldham_theme() +
    theme(legend.position='top') +
    guides(color=guide_legend(nrow=3,byrow=FALSE)) +
    guides(shape=guide_legend(nrow=2,byrow=TRUE))
dev.off()
# }}}
# nonmalignant
# {{{
me = data.frame(fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/Module_eigengenes_02-55-42.csv'))
bulkVecs = c('honeydew1', 'purple', 'plum', 'mediumorchid', 'saddlebrown', 'tan', 'lightgreen', 'lightsteelblue1', 'coral1', 'lavenderblush3', 'salmon', 'yellow', 'green', 'orange')
me = me[, c(1, which(colnames(me) %in% bulkVecs))]
me$Sample = as.numeric(gsub('.*_' , '', me$Sample))
snSlices = c(17, 53, 93, 117)
names(snSlices) = seq_along(snSlices)
inds = lapply(snSlices, function(slice) {
    order(abs(slice - me$Sample), decreasing=F)[1:4]
})
out = do.call(rbind, lapply(seq_along(inds), function(i) {
    x=inds[[i]]
    apply(me[x, ], 2, mean)
}))
# }}}
outSingle = lapply(list(9, 11, 8, 4, 12, 6), function(x) {
    temp = clust[clust==x]
    temp2 = lapply(seq_along(secs), function(k){
        sec=secs[[k]]
        (length(which(gsub('.*\\.', '', names(temp)) %in% sec))/as.numeric(secSum[k]))*100
    })
    data.frame(clone = rep(x, length(temp2)), percentage = unlist(temp2), section = c('17', '53', '93', '117'))
})
outSingle = do.call(rbind, outSingle)    
outSingle$percentage[5:8] = outSingle$percentage[5:8] + outSingle$percentage[9:12]
outP = outSingle[-seq(9,12),]
rownames(outP) = seq(1, nrow(outP))
# now combo of all nonmalignant
# {{{
test1 = reshape2::melt(out[, grep('purple$|tan$|yellow$|^green$|^orange$', colnames(out))])
test1 = data.frame(source = rep('bulk', nrow(test1)), section = test1$Var1, module = test1$Var2, percentage = test1$value*100)
test1$section = snSlices[match(test1$section, names(snSlices))]
test1$module = gsub('purple', "Module 'purple'\nAstrocyte, 1E-16", test1$module)
test1$module = gsub('green', "Module 'green'\nNeuron, 1E-36", test1$module)
test1$module = gsub('yellow', "Module 'yellow'\nMicroglia, 5E-23", test1$module)
test1$module = gsub('tan', "Module 'tan'\nOligodendrocyte, 2E-24", test1$module)
test1$module = gsub('orange', "Module 'Orange'\nEndothelial, 2E-24", test1$module)
test1$celltype = gsub('.*\n(.*?),.*', '\\1', test1$module)
outP$clone = gsub('^9$', "Clone '9'\nAstrocyte", outP$clone)
outP$clone = gsub('^11$', "Clone '11'\nOligodendrocyte", outP$clone)
outP$clone = gsub('^4$', "Clone '4'\nMicroglia", outP$clone)
outP$clone = gsub('^12$', "Clone '12'\nNeuron", outP$clone)
outP$clone = gsub('^6$', "Clone '6\nEndothelial", outP$clone)
outP$celltype = gsub('.*\n(.*?)', '\\1', outP$clone)

plot = rbind(
    data.frame(section = test1$section, frequencyVector = test1$module, frequency = test1$percentage, celltype = test1$celltype, deriv = rep('Bulk', nrow(test1))),
    data.frame(section = outP$section, frequencyVector = outP$clone, frequency = outP$percentage, celltype = outP$celltype, deriv = rep('Single nucleus', nrow(outP))))
for(ea in unique(plot$frequencyVector)){
    plot$frequency[plot$frequencyVector==ea] = scale(plot$frequency[plot$frequencyVector==ea])
}
plot$celltype = gsub('Astrocyte', 'Astrocyte, R=0.98', plot$celltype)
plot$celltype = gsub('Oligodendrocyte', 'Oligodendrocyte:1+2,\nR=0.93', plot$celltype)
plot$celltype = gsub('Microglia', 'Microglia, R=0.90', plot$celltype)
plot$celltype = gsub('Neuron', 'Neuron, R=0.99', plot$celltype)
plot$celltype = gsub('Endothelial', 'Endothelial, R=0.99', plot$celltype)
plot$section = factor(plot$section, levels=c('17', '53', '93', '117'))
plot$celltype = factor(plot$celltype, levels = c('Astrocyte, R=0.98','Oligodendrocyte:1+2,\nR=0.93','Microglia, R=0.90', 'Neuron, R=0.99', 'Endothelial, R=0.99'))
pdf(paste0(WD, 'final_nonmalig.pdf'), height=13, width=13*1.2)
ggplot(plot, aes(x =section, y=frequency, group = interaction(deriv, celltype), color=celltype, shape=deriv))+
    geom_point(size=8) +
    geom_line(size=3) +
    labs(x='Section', y='Scaled values', title='Bulk and single-nucleus non-malignant\ncell-type estimates coincide', shape = 'Source of estimate', color = 'Clone') +
    scale_color_manual(values = c('#e6ab02', '#e7298a', '#1b9e77', '#701d02','#1d4a4d')) +
    scale_y_continuous(limit = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 1.5)) +
    scale_x_discrete(expand = c(0.01,0.01)) +
    theme_classic() +
    oldham_theme() +
    theme(legend.position = 'top') +
    guides(color=guide_legend(nrow=3,byrow=FALSE)) +
    guides(shape=guide_legend(nrow=2,byrow=TRUE))
dev.off()
# }}}
# }}}
# get single-nucleus-amplicon information
# {{{
snAmp=read.csv('~/@patrick/SF10711/sn.amp.seq/malignancy.adjusted.stringent.csv')
snAmp$id = gsub('(N[0-9]*)\\.(.*)', '\\2\\.\\1', snAmp$id)
snAmp=snAmp[match(colnames(exprD), snAmp$id),]
snAmp$malignant = gsub('0', 'Nonmalignant', snAmp$malignant)
snAmp$malignant = gsub('1', 'Malignant', snAmp$malignant)
snAmp$malignant = gsub('UNK', 'No data', snAmp$malignant)
# }}}
# panel E - copykat data
# {{{
copyDat = read.csv('/home/patrick/@patrick/SF10711/sn.figs/copykat_cytoband_summary.csv')[,-seq(46,49)]
rownames(copyDat) = copyDat[,1]
copyDat = copyDat[,-1]
copyDat = apply(copyDat, 2, as.numeric)
copyDat = as.matrix(copyDat)
copyD = hclust(dist(copyDat[,c(8,9, 14, 16, 36)], method = 'euclidean'), method = 'ward.D')
copyDClust = cutree(copyD, k=4)
copyDClust = c('No CNV detected', 'CNV cluster 1', 'CNV cluster 2', 'CNV cluster 3')[match(copyDClust, c(2, 3, 1, 4))]
colnames(copyDat) = gsub('X', '', colnames(copyDat))
copyDat[copyDat<0.05 & copyDat>(-0.05)]=0
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
clustF = names(clustConv)[match(clust, clustConv)]
clustF = factor(clustF, levels = c('Clone 1', 'Clone 2', 'Clone 3', 'Clone 4:1', 'Clone 4:2', 'Clone 5', 'Astrocytes', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons', 'Endothelial cells'))
clustCol = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a', '#90af41', '#a65729', '#e6ab02', '#e7298a', '#e729e6', '#1b9e77', '#701d02', '#1d4a4d')
names(clustCol) = c('Clone 1', 'Clone 2', 'Clone 3', 'Clone 4:1', 'Clone 4:2', 'Clone 5', 'Astrocytes', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons', 'Endothelial cells')
rha = rowAnnotation(
    'Amplicon-seq' = snAmp$malignant,
    'snRNA-seq clusters' = clustF,
     col = list('Amplicon-seq' = c('Nonmalignant' = '#2166ac', 'Malignant' = '#b2182b', 'No data' = '#999999'),
                'snRNA-seq clusters' = clustCol)
)
copyDatCounts = apply(copyDat, 1, function(x) sum(abs(x), na.rm=T))
row_ha = rowAnnotation('Total calls' = anno_barplot(copyDatCounts))
pdf(paste0(WD, 'copykat_figure.pdf'), height=7, width = 14)
Heatmap(copyDat,
    name = 'Log ratio',
    cluster_columns = FALSE,
    column_names_side = "top",
    column_names_centered = TRUE,
    column_names_gp = gpar(fontfamily = 'NimbusSan', cex = 0.8),
    row_title_gp = gpar(fontfamily = 'NimbusSan', cex = 1.5),
    row_dend_width = unit(2, 'cm'),
    cluster_rows = copyD,
    row_title = NULL,
    column_names_rot = 0, 
    row_split = 4,
    left_annotation = rha,
    right_annotation = row_ha,
    col = colorRamp2(breaks = c(-.1, 0, .1),
                            colors = c('#2166ac', '#f7f7f7', '#b2182b')
                            )
)
dev.off()
# }}}
# panel C - heatmap showing the top differentially expressed genes for each clone. On the left hand side, will show colored vectors for copykat and amplicon calls, total expression levels from the cell, as well as the dendrogram. Maker labels will be on the top.
# {{{
# Identify top differentially expressed genes based on t-value for each clone
# {{{
## for this one do wilcoxn rank sum test
out = future_lapply(unique(clust), function(clu){
    ind = (clu == clust)
    future_apply(exprD, 1, function(x){
        wilcox.test(x[ind], x[!ind], alternative='greater')$p.value
    })
})
exprNVec = c(2, 10, 9, 7, 1,6,  5, 3, 11, 8, 4, 12)
exprNVec = exprNVec[seq(length(exprNVec), 1)]
exprN = future_lapply(exprNVec, function(i){
    ind = order(out[[i]], decreasing=F)[1:15000]
    try = data.frame(exprGene[ind], exprD[ind,])
    try = try[!(is.na(try[,1])),]
    if(length(grep('-AS|LINC', try[,1])>0)){
        try = try[-grep('-AS|LINC', try[,1]), ]
    }
    try[1:15000,]
})
# manually selected genes that are protein coding (excluding non-coding genes)
# {{{
exprN[[1]] = exprN[[1]][c(1,3, 4,5,6),]
exprN[[2]] = exprN[[2]][1:5,]
exprN[[3]] = exprN[[3]][c(3,4,5,6,7),]
exprN[[4]] = exprN[[4]][1:5,]
exprN[[5]] = exprN[[5]][1:5,]
exprN[[6]] = exprN[[6]][c(2,3,5, 6, 8),]
exprN[[7]] = exprN[[7]][c(1,2,3,4,7),]
exprN[[8]] = exprN[[8]][c(1,2,3,5, 18),]
exprN[[9]] = exprN[[9]][c(1, 3, 4,5,6),]
exprN[[10]] = exprN[[10]][c(1,2,3,4,18),]
exprN[[11]] = exprN[[11]][2:6,]
exprN[[12]] = exprN[[12]][1:5,]
# exprN = lapply(exprN, function(x) x[1:5,])
# }}}
exprN = do.call(rbind, exprN)
plotG = exprN[,1]
exprN = exprN[,-1]
plot = t(as.matrix(log2(exprN+1)))
plot = apply(plot, 2, scale)
colnames(plot) = plotG
# }}}
# Get total umi counts for each nucleus
# {{{
counts = apply(exprC[, -c(1,2)], 2, sum)
# }}}
# Get copkykat cnv information
# {{{
# Add vector for copykat on side (from: /@patrick/SF10711/cnv.analysis/sn.rna.seq/copykat/21_03_03/21_03_09_SF10711_hi_sens_nomag_base_copykat/log_21_03_03.R)
cna =fread('~/@patrick/SF10711/cnv.analysis/sn.rna.seq/copykat/base_copykat_CNA_results.txt')
cna_t_vals=melt(cna[,-c(1,2,3)])
chr_select=seq(1,22)
cna_t=cna[which(cna$chrom %in% chr_select),]
cna_clust=cna[which(cna$chrom %in% c(1,2,7,8,9,10)),-c(1,2,3)]
cna_clust=apply(cna_clust, 2, as.numeric)
cna_clust[intersect(which(cna_clust<mean(cna_t_vals$value)+1*sd(cna_t_vals$value)), which(cna_clust>as.numeric(mean(cna_t_vals$value)-1*sd(cna_t_vals$value))))]=0
cna_hclust=hclust(parDist(t(cna_clust), threads=15, method='euclidean'),method='ward.D2')
funCNA=function(cna_clust2, knu){list(cluster=cutree(hclust(parDist(t(cna_clust2), threads=15, method='euclidean'),method='ward.D2'), knu))}
clusters=cutree(cna_hclust, k=3)
names(clusters) = gsub('(N[0-9]*)\\.(.*)', '\\2\\.\\1', names(clusters))
clusters = clusters[match(colnames(exprD), names(clusters))]
clusters = gsub('^2$', 'No CNV detected', clusters)
clusters = gsub('^1$', 'CNV cluster 1', clusters)
clusters = gsub('^3$', 'CNV cluster 2', clusters)
clusters = as.factor(clusters)
# }}}
# Infer CNV because copykat data is discrepeant with snamp-seq, will also include data from infercnv
# {{{
inferDend = read.dendrogram('~/@patrick/SF10711/cnv.analysis/sn.rna.seq/infercnv/infercnv.preliminary.observations_dendrogram.txt')
inferClust = cutree(inferDend, k=3)
indeces=c("TAAGGCGA", "CGTACTAG", "AGGCAGAA", "TCCTGAGC", "GGACTCCT", "TAGGCATG", "CTCTCTAC", "CAGAGAGG", "GCTACGCT", "CGAGGCTG", "AAGAGGCA", "GTAGAGGA")
names(indeces)=paste0('N', seq(701,712))
for(ea in indeces){
    i = which(indeces == ea)
    convInd = which(ea == substr(names(inferClust), 1, 8))
    temp = paste('X', substr(names(inferClust)[convInd], 9,20), names(indeces)[i], sep='.')
    names(inferClust)[convInd] = temp
}
names(inferClust) = gsub('^X\\.', '', names(inferClust))
inferClust = inferClust[match(snAmp$id, names(inferClust))]
inferClust = gsub('^3$', 'No CNV detected ', inferClust)
inferClust = gsub('^2$', 'CNV cluster 1 ', inferClust)
inferClust = gsub('^1$', 'CNV cluster 2 ', inferClust)
# }}}
# get CaSpER clusters
# {{{
casperClust = read.csv('~/@patrick/SF10711/cnv.analysis/sn.rna.seq/casper/CaSpER_clusters.csv')[,2]
casperClust = gsub('2', 'No CNV detected  ', casperClust)
casperClust = gsub('1', 'CNV detected', casperClust)
# }}}
# get mitosis signature
# {{{
library('flashClust')
source("/home/patrick/code/git/GSEA_generic/GSEAfxsV3.R")
broadSets=getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
geneSelec = geneIds(broadSets[18950])
geneInds = match(unlist(geneSelec), exprGene)
geneInds = geneInds[!(is.na(geneInds))]
mitosisVec = as.numeric(apply(exprD[geneInds,], 2, mean))
# }}} 
# Create heatmap (complex heatmap)
# {{{
row_ha = rowAnnotation('Unique reads' = anno_barplot(counts), 'Mitosis signature' = anno_barplot(mitosisVec))
row_ha_malig = rowAnnotation(
    'Amplicon-seq' = snAmp$malignant,
    'CopyKat CNV' = copyDClust,
    'InferCNV' = inferClust, 
    'CaSpER' = casperClust,
    col = list('Amplicon-seq' = c('Nonmalignant' = '#2166ac', 'Malignant' = '#b2182b', 'No data' = '#999999'), 
              'CopyKat CNV' = c('No CNV detected' = '#2166ac', 'CNV cluster 1' = '#b2182b', 'CNV cluster 2' = '#4daf4a', 'CNV cluster 3' = '#ffff33'),
              'InferCNV' = c('No CNV detected ' = '#2166ac', 'CNV cluster 1 ' = '#b2182b', 'CNV cluster 2 ' = '#4daf4a'),
              'CaSpER' = c('No CNV detected  ' = '#2166ac', 'CNV detected' = '#b2182b')
              )
)
dendCol = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a', '#90af41', '#a65729', '#e6ab02', '#e7298a', '#e729e6', '#1b9e77', '#701d02', '#1d4a4d')[c(11, 10, 9, 8, 4, 5, 12, 2, 3, 7, 6, 1)]
dendP = as.dendrogram(distHclust)
dendP = color_branches(dendP, 
                        k = 12,
                        col = dendCol)
pdf(paste0(WD, 'sn_heatmap.pdf'), height=7, width = 16)
Heatmap(plot,
    name = 'Scaled log2\nexpression',
    cluster_rows = dendP,
    col = colorRamp2(breaks = c(-4, 0, 4),
                            colors = c('#2166ac', '#f7f7f7', '#b2182b')
                            ),
    row_split = 12,
    row_dend_width = unit(2, 'cm'),
    row_title = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Endothelial cells', 'Clone 4:2', 'Clone 4:1','Oligodendrocytes:2', 'Oligodendrocyte:1', 'Microglia', 'Neurons')[seq(12,1)],
    row_title_gp = gpar(col = dendCol, fontfamily = 'NimbusSan', cex = 1.5),
    row_title_rot = 0,
    column_title = NULL,
    cluster_columns = FALSE,
    column_names_rot = 45, 
    column_names_side = "top",
    column_names_gp = gpar(col = c(rep(dendCol, ea =5), 'red'), fontfamily = 'NimbusSan', cex = 0.8),
    column_split = c(rep(seq(1,12), ea = 5), 13),
    right_annotation = row_ha,
    left_annotation = row_ha_malig
)
dev.off()
# }}}
# }}}
# panel D1 - umap with clonal definitions
# {{{
exprDbak =exprD
exprD = log2(exprD+1)
exprPR = prcomp(t(exprD))
set.seed(15)
exprUmap = uwot::umap(exprPR$x[,1:30], 
                        n_neighbors = 30, 
                        spread = 3, 
                        min_dist = 2, 
                        n_threads = 1, 
                        n_sgd_threads = 1, 
                        metric = 'correlation', 
                        n_epochs = 1000)
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocytes', 'Clone 3', 'Clone 2', 'Clone 4:2', 'Endothelial cells', 'Clone 4:1', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons')
clustF = names(clustConv)[match(clust, clustConv)]
dendCol = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a', '#90af41', '#a65729', '#e6ab02', '#e7298a', '#e729e6', '#1b9e77', '#701d02', '#1d4a4d')# [c(11, 10, 9, 8, 4, 5, 12, 2, 3, 7, 6, 1)]
colnames(exprUmap) = c('umap_1', 'umap_2')
tempDt = data.table(exprUmap, col = factor(clustF, levels = c('Clone 1', 'Clone 2', 'Clone 3', 'Clone 4:1', 'Clone 4:2', 'Clone 5', 'Astrocytes', 'Oligodendrocytes:1', 'Oligodendrocytes:2', 'Microglia', 'Neurons', 'Endothelial cells')))
pdf(paste0(WD, 'umap_hclust.pdf'), height=13, width=13*1.2)
ggplot(tempDt, aes(x = umap_1, y = umap_2, color = col)) +
    geom_point(size = 5) +
    scale_color_manual(values = dendCol)+
    theme_classic() +
    labs(title = 'UMAP plot of snRNA-seq clusters',
         x = 'UMAP dim 1',
         y = 'UMAP dim 2',
         color = 'Clusters') +
    oldham_theme() +
    theme(legend.position = 'top') +
    guides(color=guide_legend(nrow=4,byrow=TRUE))
dev.off()
# }}}
# panel D2 - umap with mutation defintions
# {{{
clustF = snAmp$malignant[match(rownames(exprUmap), snAmp$id)]
clustF = snAmp$malignant[match(rownames(exprUmap), snAmp$id)]
tempDt = data.table(exprUmap, 
            col = factor(clustF, levels = c('Malignant', 'Nonmalignant', 'No data')))
dendCol = c('#b2182b', '#2166ac', '#878787')
pdf(paste0(WD, 'umap_hclust_malig.pdf'), height=13, width=13*1.2)
ggplot(tempDt, aes(x = umap_1, y = umap_2, color = col)) +
    geom_point(size = 5) +
    scale_color_manual(values = dendCol)+
    theme_classic() +
    labs(title = 'UMAP plot of snRNA-seq clusters',
         x = 'UMAP dim 1',
         y = 'UMAP dim 2',
         color = 'Clusters') +
    oldham_theme() +
    theme(legend.position = 'top') +
    guides(color=guide_legend(nrow=1,byrow=TRUE))
dev.off()
# }}}
# panel D3 - umap with section definitions
# {{{
secs = list(c('N701','N702', 'N709'), c('N712', 'N703', 'N704'), c('N705', 'N706', 'N710'), c('N707', 'N708', 'N711'))
names(secs) = c('s17', 's53', 's93', 's117')
clustF = rep(NA, nrow(exprUmap))
for(i in seq_along(secs)){
        x = secs[[i]]
        clustF[grep(paste0('.*\\.(', paste0(x, collapse = '|'), ')$', collapse = ''), rownames(exprUmap))] = names(secs)[i]
}
tempDt = data.table(exprUmap,
                    col = factor(clustF, levels = c('s17', 's53', 's93', 's117')))
out$section = gsub('s14', 's17', out$section)
tempDt = data.table(exprUmap, col = factor(out$section[match(rownames(exprUmap), out$nuclei)], levels = c('s17', 's53', 's93', 's117')))

dendCol = c('#bdc9e1', '#67a9cf', '#1c9099', '#016c59')
pdf(paste0(WD, 'umap_hclust_section.pdf'), height=13, width=13*1.2)
ggplot(tempDt, aes(x = umap_1, y = umap_2, color = col)) +
    geom_point(size = 5) +
    scale_color_manual(values = dendCol)+
    theme_classic() +
    labs(title = 'UMAP plot of snRNA-seq clusters',
         x = 'UMAP dim 1',
         y = 'UMAP dim 2',
         color = 'Clusters') +
    oldham_theme() +
    theme(legend.position = 'top') +
    guides(color=guide_legend(nrow=1,byrow=TRUE))
dev.off()
# }}}
# panel E - copykat data
# {{{
copyDat = read.csv('/home/patrick/code/git/SF10711_single_nucleus/copykat_cytoband_summary.csv')[,-seq(46,49)]
rownames(copyDat) = copyDat[,1]
copyDat = copyDat[,-1]
copyDat = apply(copyDat, 2, as.numeric)
copyDat = as.matrix(copyDat)
copyD = hclust(dist(copyDat, method = 'euclidean'), method = 'ward.D2')
colnames(copyDat) = gsub('X', '', colnames(copyDat))

clustConv = c(2, 9, 8, 6, 1, 3, 5, 7, 4, 10)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocyte', 'Clone 3', 'Unknown', 'Clone 4', 'Clone 4','Oligodendrocyte', 'Microglia', 'Neuron')
clustF = names(clustConv)[match(clust, clustConv)]
clustF = factor(clustF, levels = c('Clone 1', 'Unknown', 'Clone 3', 'Clone 4', 'Clone 5', 'Astrocyte', 'Oligodendrocyte', 'Microglia', 'Neuron'))
clustCol = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a','#a65729', '#e6ab02', '#e7298a', '#1b9e77', '#701d02')
names(clustCol) = c('Clone 1', 'Unknown', 'Clone 3', 'Clone 4', 'Clone 5', 'Astrocyte', 'Oligodendrocyte', 'Microglia', 'Neuron')
rha = rowAnnotation(
    'Amplicon-seq' = snAmp$malignant,
    'snRNA-seq clusters' = clustF,
     col = list('Amplicon-seq' = c('Nonmalignant' = '#2166ac', 'Malignant' = '#b2182b', 'No data' = '#999999'),
                'snRNA-seq clusters' = clustCol)
)
copyDatCounts = apply(copyDat, 1, function(x) sum(abs(x)))
row_ha = rowAnnotation('Total calls' = anno_barplot(copyDatCounts))
pdf('copykat_figure.pdf', height=7, width = 14)
Heatmap(copyDat,
    name = 'Log ratio',
    cluster_columns = FALSE,
    column_names_side = "top",
    column_names_centered = TRUE,
    column_names_gp = gpar(fontfamily = 'NimbusSan', cex = 0.8),
    row_title_gp = gpar(fontfamily = 'NimbusSan', cex = 1.5),
    row_dend_width = unit(2, 'cm'),
    cluster_rows = copyD,
    row_title = NULL,
    column_names_rot = 0, 
    row_split = 5,
    left_annotation = rha,
    right_annotation = row_ha,
    col = colorRamp2(breaks = c(-.15, 0, .15),
                            colors = c('#2166ac', '#f7f7f7', '#b2182b')
                            )
)
dev.off()
# }}}
# panel F - trajectory
# {{{
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocyte', 'Clone 3', 'Clone 2', 'Clone 4:1', 'Endothelial\ncells', 'Clone 4:2', 'Oligodendrocyte:1', 'Oligodendrocyte:2', 'Microglia', 'Neuron')
nonmalignant = c(9, 11, 8, 4, 12, 6)
maligVec = names(clust)[!(clust %in% nonmalignant)]
clustMalig = clust[!(clust %in% nonmalignant)]

exprDTraj = log2(exprD[,which(colnames(exprD) %in% maligVec)]+1)
trajPC = prcomp(t(exprDTraj))$x[,1:15]
set.seed(15)
trajUMAP = uwot::umap(trajPC, n_neighbors=20, spread=2, min_dist=1, n_threads=1, n_sgd_threads=1, metric='correlation')
trajUMAP = uwot::umap(trajPC, n_neighbors=40, spread=5, min_dist=1, n_threads=1, n_sgd_threads=1, metric='correlation')
lin1 <- getLineages(trajUMAP, clustMalig, start.clus = NULL, dist.method = 'simple')
crv1 <- getCurves(lin1)

clustF = names(clustConv)[match(clustMalig, clustConv)]

dendCol = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a', '#90af41', '#a65729', '#e6ab02', '#e7298a', '#e729e6', '#1b9e77', '#701d02', '#1d4a4d')# [c(11, 10, 9, 8, 4, 5, 12, 2, 3, 7, 6, 1)]
colnames(trajUMAP) = c('umap_1', 'umap_2')
trajUMAPP = data.table(trajUMAP, col = factor(clustF, levels = c('Clone 1', 'Clone 2', 'Clone 3', 'Clone 4:1', 'Clone 4:2', 'Clone 5')))
write.csv(trajUMAPP, file = paste0(WD,'umap_for_trajectory_analysis.csv'))

nodes = t(as.data.frame(igraph::get.vertex.attribute(lin1@metadata$mst, 'coordinates')))
# igraph::get.edgelist(lin1@metadata$mst)
rownames(nodes) = igraph::get.vertex.attribute(lin1@metadata$mst, 'name')
colnames(nodes) = c('dim1', 'dim2')

pdf(paste0(WD, 'trajectory.pdf'),height=13, width=17)
ggplot(trajUMAPP, aes(x = umap_1, y = umap_2, color = col)) +
    geom_point(size = 5) +
    scale_color_manual(values = dendCol)+
    geom_segment(aes(x = nodes[2,1], y = nodes[2,2], xend = nodes[1,1], yend = nodes[1,2]), 
                    lineend = 'round', 
                    color = 'black', arrow = arrow(length = unit(0.3, "inches")),
                    size = 3) +
    geom_segment(aes(x = nodes[1,1], y = nodes[1,2], xend = nodes[5,1], yend = nodes[5,2]), 
                    lineend = 'round',arrow = arrow(length = unit(0.3, "inches")), 
                    color = 'black', 
                    size = 3) +
    geom_segment(aes(x = nodes[2,1], y = nodes[2,2], xend = nodes[3,1], yend = nodes[3,2]), 
                    lineend = 'round', arrow = arrow(length = unit(0.3, "inches")),
                    color = 'black', 
                    size = 3) +
    geom_segment(aes(x = nodes[3,1], y = nodes[3,2], xend = nodes[4,1], yend = nodes[4,2]), 
                    lineend = 'round', arrow = arrow(length = unit(0.3, "inches")),
                    color = 'black', 
                    size = 3) +
    geom_segment(aes(x = nodes[2,1], y = nodes[2,2], xend = nodes[6,1], yend = nodes[6,2]), 
                    lineend = 'round', arrow = arrow(length = unit(0.3, "inches")),
                    color = 'black', 
                    size = 3) +
    labs(title = 'UMAP plot of malignant cells with trajectory analysis',
         x = 'UMAP dim. 1',
         y = 'UMAP dim 2',
         color = 'Clusters') +
    theme_classic() +
    oldham_theme()+
    theme(legend.position='top')
dev.off()

nonmalignant = c(8, 7, 4, 10)
controlSampleIDs = clusters$X[clusters$x %in% nonmalignant]
maligVec = clusters$x[!(clusters$x %in% nonmalignant)]
maligSampleIDs = clusters$X[!(clusters$x %in% nonmalignant)]

trajUMAP = read.csv('~/@patrick/SF10711/sn.figs/umap_for_trajectory_analysis.csv')
clustConv = c(2, 10, 9, 7, 1, 5, 6, 3, 11, 8, 4, 12)
names(clustConv) = = c('Clone 1', 'Clone 5', 'Astrocyte', 'Clone 3', 'Clone 2', 'Clone 4:1', 'Endothelial\ncells', 'Clone 4':2, 'Oligodendrocyte:1', 'Oligodendrocyte:2', 'Microglia', 'Neuron')
dendColMalig = dendCol[c(1,2,3,4,5,6)]
clustersMalig = clusters$x[!(clusters$x %in% c(9, 11, 8, 4, 12, 6))]
 trajUMAP$col = names(clustConv)[match(clustersMalig, clustConv)]
trajUMAP$col = factor(trajUMAP$col, levels = c('Clone 1', 'Unknown', 'Clone 3', 'Clone 4 : 1', 'Clone 4 : 2', 'Clone 5'))
colnames(trajUMAP) = c('umap_1', 'umap_2', 'col')

trajUMAP = umap(trajPC, n_neighbors=20, spread=2, min_dist=1, n_threads=15, n_sgd_threads=15, metric='correlation')
clustConv = c(2, 9, 8, 6, 1, 3, 5, 7, 4, 10)
names(clustConv) = c('Clone 1', 'Clone 5', 'Astrocyte', 'Clone 3', 'Unknown', 'Clone 4 : 1', 'Clone 4 : 2','Oligodendrocyte', 'Microglia', 'Neuron')
clustF = names(clustConv)[match(maligVec, clustConv)]
dendCol = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a','#a65729', '#e6ab02', '#e7298a', '#1b9e77', '#701d02')
colnames(trajUMAP) = c('umap_1', 'umap_2')
trajUMAPP = data.table(trajUMAP, col = factor(clustF, levels = c('Clone 1', 'Unknown', 'Clone 3', 'Clone 4 : 1', 'Clone 4 : 2', 'Clone 5', 'Astrocyte', 'Oligodendrocyte', 'Microglia', 'Neuron')))
trajUMAP= trajUMAPP
lin1 <- getLineages(trajUMAP[,c(1,2)], trajUMAP$col, start.clus = 'Clone 1', dist.method = 'simple')
nodes = t(as.data.frame(igraph::get.vertex.attribute(lin1@metadata$mst, 'coordinates')))
igraph::get.edgelist(lin1@metadata$mst)
rownames(nodes) = igraph::get.vertex.attribute(lin1@metadata$mst, 'name')
colnames(nodes) = c('dim1', 'dim2')
dendCol = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a', '#2b6129', '#a65729', '#e6ab02', '#e7298a', '#1b9e77', '#701d02')
pdf(paste0(WD, 'trajectory.pdf'),height=13, width=13*1.1)
ggplot(trajUMAP, aes(x = umap_1, y = umap_2, color = col)) +
    geom_point(size = 5) +
    scale_color_manual(values = dendCol)+
    geom_segment(aes(x = nodes[2,1], y = nodes[2,2], xend = nodes[6,1], yend = nodes[6,2]), # clone 1 to clone 5
                    lineend = 'round', 
                    color = 'black', arrow = arrow(length = unit(0.3, "inches")),
                    size = 3) +
    geom_segment(aes(x = nodes[2,1], y = nodes[2,2], xend = nodes[1,1], yend = nodes[1,2]), # clone 1 to clone 2
                    lineend = 'round',arrow = arrow(length = unit(0.3, "inches")), 
                    color = 'black', 
                    size = 3) +
    geom_segment(aes(x = nodes[1,1], y = nodes[1,2], xend = nodes[5,1], yend = nodes[5,2]), # clone 2 to clone 3
                    lineend = 'round', arrow = arrow(length = unit(0.3, "inches")),
                    color = 'black', 
                    size = 3) +
    geom_segment(aes(x = nodes[2,1], y = nodes[2,2], xend = nodes[3,1], yend = nodes[3,2]), # clone 2 to clone 3
                    lineend = 'round', arrow = arrow(length = unit(0.3, "inches")),
                    color = 'black', 
                    size = 3) +
    geom_segment(aes(x = nodes[3,1], y = nodes[3,2], xend = nodes[4,1], yend = nodes[4,2]), #clone 4:1 to clone 4:2
                    lineend = 'round', arrow = arrow(length = unit(0.3, "inches")),
                    color = 'black', 
                    size = 3) +
    labs(title = 'UMAP plot of malignant cells\nwith trajectory analysis',
         x = 'UMAP dim 1',
         y = 'UMAP dim 2',
         color = 'Clusters') +
    theme_classic() +
    oldham_theme()+ 
    theme(legend.position = 'top')
dev.off()
# }}}
# enrichments of celltypes malignant and nonmalignant
# calculate  sensitivity and specificity
# sensitivity definition: (true positives) / (true positives + false negatives)
# {{{
sensitivity = function(truth, test){
   temp = data.frame(truth,test)
   tp = length(which(apply(temp, 1, function(x) x[1] == 1 & x[2] == 1)))
   fn = length(which(apply(temp, 1, function(x) x[1] == 1 & x[2] == 0)))
   return(tp/(tp+fn))
}
# specificity definition: (true negative) / (true negatives + false positives)
specificity = function(truth, test){
   temp = data.frame(truth,test)
   tn = length(which(apply(temp, 1, function(x) x[1] == 0 & x[2] == 0)))
   fp = length(which(apply(temp, 1, function(x) x[1] == 0 & x[2] == 1)))
   return(tn/(tn+fp))
}
# accuracy definition: (TP + TN) / (TP + FP + TN + FN)
accuracy = function(truth, test){
   temp = data.frame(truth,test)
   tp = length(which(apply(temp, 1, function(x) x[1] == 1 & x[2] == 1)))
   fp = length(which(apply(temp, 1, function(x) x[1] == 0 & x[2] == 1)))
   tn = length(which(apply(temp, 1, function(x) x[1] == 0 & x[2] == 0)))
   fn = length(which(apply(temp, 1, function(x) x[1] == 1 & x[2] == 0)))
   return((tp + tn)/(tp + tp + tn + fn))
}

# snAmp-seq ground truth
truth = snAmp$malignant
truth = gsub('No data', NA, truth)
truth = truth == 'Malignant'
# CopyKat
sensitivity(truth, copyDClust != 'No CNV detected')
specificity(truth, copyDClust != 'No CNV detected')
accuracy(truth, copyDClust != 'No CNV detected')
# InferCNV
sensitivity(truth, inferClust != 'No CNV detected ')
specificity(truth, inferClust != 'No CNV detected ')
accuracy(truth, inferClust != 'No CNV detected ')
# CaSpER
sensitivity(truth, casperClust != 'No CNV detected  ')
specificity(truth, casperClust != 'No CNV detected  ')
accuracy(truth, casperClust != 'No CNV detected  ')
# }}}
# get FM clone genesets and nonmalignant genesets
# {{{
# we can instead get the genes from the dataset
clones = c('Endothelial_orange', 'Endothelial_plum1', 'Endothelial_mediumorchid', 'Endothelial_brown4', 'Clone4_ivory')
# orange is endothelial cells , number 2, *supplement*
# ivory is cluster 4, which goes to 1 and 3
bulkGenes = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
clones=c('orange', 'ivory')
bulkGenes = lapply(clones, function(clone) bulkGenes$Gene[which(bulkGenes$'TopModPosFDR_0.12' == clone)])
names(bulkGenes) = c('Orange - endothelial cell', 'Ivory - clone 4') 
# }}}
# get cluster definitions from hierarchical clustering
# {{{
attr(distDf, 'Labels') =paste0('x', seq(1, length(attr(distDf, 'Labels'))))
distHclust = hclust(distDf, method = 'ward.D')
dendP = as.dendrogram(distHclust)
clust = cutree(distHclust, k=12)
names(clust) = colnames(expr)[-c(1,2)]
# }}}
# panel G, Then reuse code to get distributions of t-values
# {{{
tValues = function(exprMat, inds){
    apply(exprMat, 1, function(i) t.test(as.numeric(i[inds]), as.numeric(i[!inds]), paired=F)$statistic)
}
# names(set1) = c('Red - clone 1', 'Black - clone 3', 'Ivory - clone 4', 'Lightcyan - clone 5', 'Honeydew1 - astrocyte', 'Tan - oligodendrocyte', 'Coral1 - microglia', 'Lavenderblush3 - neuron')
tVal = future_lapply(bulkGenes, function(set){
    ind2 = exprGene %in% set
    future_lapply(c(3,5,6), function(clu){
        ind = (clust == clu)
        print(clu)
        list(tValues(exprD[ind2,], inds = ind), tValues(exprD[!ind2,], inds = ind))
    })
})
sampN = 1000
adOut = future_lapply(tVal, function(y){
    future_lapply(y, function(x){
        wilcox.test(x[[1]], x[[2]][sample(seq(1, length(x[[2]])), 100)], exact = F, alternative = 'greater', correct = F)$p.value*10*8
    })
})
adOut = do.call(rbind , adOut)
colnames(adOut) = paste('Cluster', c('4 : 1', '4 : 2 : 1', '4 : 2 : 2'))

# plot
# {{{
sigVec = list(3, c(1,2))
annoVec = lapply(sigVec, function(x) rep('***', length(x)))
colTitle = c('#ff7f00', '#e6ab02', '#dd95e8', '#e7298a', '#4daf4a', '#1b9e77', '#a65729', '#701d02')
rowN=seq(1,3)
names(rowN) = c('4 : 1', '4 : 2 : 1', '4 : 2 : 2')
plotZ = lapply(seq(1,2), function(i){
    h = i
    try = lapply(seq_along(tVal[[1]]), function(j){
        plot = tVal[[i]][[j]]
        plot2 = reshape2::melt(plot)
        plot2$L1 = gsub(1, 'Gene set', plot2$L1)
        plot2$L1 = gsub(2, 'All other genes', plot2$L1)
        plot2$Clone = j
        colnames(plot2) = c('T_value', 'Set', 'Clone')
        plot2
    })
    try = do.call(rbind, try)
    try$Set = factor(try$Set, levels = c('Gene set', 'All other genes'))
    try$Clone= factor(try$Clone, levels = seq(1,10))
    try$T_value[try$T_value>10] = NA
    try$T_value[try$T_value<(-10)] = NA
    try$Clone = names(rowN)[match(try$Clone, rowN)]
    ggplot(try, aes(x = Clone, y = T_value, fill = Set))+
        geom_violin(scale='area') +
        scale_y_continuous(limits=c(-15, 17), breaks = seq(-10, 10, 10))+
        labs(y = '', x = '', title = names(tVal)[i], fill='') +
        geom_hline(yintercept=0, color='grey50', size=0.6)+
        scale_fill_manual(values=c('white', 'black')) +
        geom_signif(annotations = annoVec[[h]],
            xmin = sigVec[[h]]-.3, xmax = sigVec[[h]]+.3,  y_position = 11, tip_length = 0, textsize = 15, vjust=.5, size = 1, family = 'NimbusSan'
        ) +
        theme_classic() +
        oldham_theme() +
        theme(axis.text.x = element_text(size=40, color='black',  family='NimbusSan', margin=margin(t=10)),
                axis.text.y = element_text(size=40, color='black',  family='NimbusSan', margin=margin(r = 15)),
                panel.grid.major.x = element_line(colour = "grey50", size = 0.6),
                plot.title = element_text(colour = colTitle[h]),
                plot.margin = unit(c(0,0,0,0), "cm"), 
                legend.position = 'none')
})
pdf(paste(WD, 'clone3_5_supp.pdf'), height =, 13*1, width = 13)
plotL=list(grobs=plotZ, ncol=1, 
            bottom = textGrob('Clusters',  gp = gpar(fontsize = 65, fontface='bold', fontfamily='NimbusSan')),
            left = textGrob('T-values',  gp = gpar(fontsize = 65, fontface='bold', fontfamily='NimbusSan'), rot = 90)
)
margin = theme()
grid.arrange(do.call(arrangeGrob, plotL), padding = unit(c(0,0,0,0), "cm"))
dev.off()
# }}}
# }}}
# panel H and I, relative abuandance for each of the clones
# {{{
# nonmalignant single cell
# {{{
secs = list(c('N701','N702', 'N709'), c('N712', 'N703', 'N704'), c('N705', 'N706', 'N710'), c('N707', 'N708', 'N711'))
secSum = lapply(secs, function(sec) length(grep(paste(sec, collapse = '|'), names(clust))))
clonesInd = c(3,5,6)
outP = lapply(seq_along(clonesInd), function(j){
    i = clonesInd[j]
    temp = clust[clust==i]
    temp2 = lapply(seq_along(secs), function(k){
        sec=secs[[k]]
        (length(which(gsub('.*\\.', '', names(temp)) %in% sec))/as.numeric(secSum[k]))*100
    })
    data.frame(clone = rep(i, length(temp2)), percentage = unlist(temp2), section = c('17', '53', '93', '117'))
})
outP = do.call(rbind, outP)
outP = data.frame(source =rep('single_cell', nrow(outP)), outP)
outP$clone = gsub('^3$', 'clone.4:1', outP$clone)
outP$clone = gsub('^5$', 'clone.4:2:1', outP$clone)
outP$clone = gsub('^6$', 'clone.4:2:2', outP$clone)
# }}}
# malignant bulk
# {{{
freq = read.csv('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv')
freq$ampseq.sample = as.numeric(gsub('.*_' , '', freq$ampseq.sample))
snSlices = c(17, 53, 93, 117)
inds = lapply(snSlices, function(slice) {
    order(abs(slice - freq$ampseq.sample), decreasing=F)[1:4]
})
weights = c(
    length(grep("N701|N702|N709", colnames(expr))),
    length(grep("N712|N703|N704", colnames(expr))),
    length(grep("N705|N706|N710", colnames(expr))),
    length(grep("N707|N708|N711", colnames(expr)))
)
out = do.call(rbind, lapply(seq_along(inds), function(i) {
    x=inds[[i]]
    apply(freq[x, ], 2, mean)
}))
out = rbind(out, apply(out, 2, function(x) {
    mean(
        c(rep(x[1], weights[1]),
        rep(x[2], weights[2]),
        rep(x[3], weights[3]),
        rep(x[4], weights[4]))
        )
    })
)
out = reshape2::melt(out[-5,c(1, 6)])
out$Var1 = gsub('^1$','17', out$Var1)
out$Var1 = gsub('^2$','53', out$Var1)
out$Var1 = gsub('^3$','93', out$Var1)
out$Var1 = gsub('^4$','117', out$Var1)
out = data.frame(source=rep('Bulk', nrow(out)),clone = out$Var2, percentage=out$value*100, section=out$Var1)[5:8,]
maligBulk = out
# }}}
# panel I, plot malignant 
# {{{
outMalig = rbind(outP[seq(1,8),], maligBulk)
rownames(outMalig)=seq_len(nrow(outMalig))
cor(outMalig$percentage[1:4],outMalig$percentage[9:12])
cor(outMalig$percentage[1:4] + outMalig$percentage[5:8],outMalig$percentage[9:12])
outMalig = rbind(outMalig, data.frame(source = rep('Single-nucleus', 4), 
                                                    clone = rep('Single-nucleus clone\n4:1 + 4:2:1, R=0.98', 4), 
                                                    percentage = outMalig$percentage[1:4] + outMalig$percentage[5:8],
                                                    section = c(17, 53, 93, 117)
                            ))
cor(outMalig$percentage[13:16],outMalig$percentage[9:12])
outMalig$source = gsub('single_cell', 'Single-nucleus', outMalig$source)
outMalig$section = factor(outMalig$section, levels = c(17, 53, 93, 117))
outMalig$clone = gsub('clone.4:1', 'Single-nucleus clone\n4:1, R=0.85', outMalig$clone, fixed = T)
outMalig$clone = gsub('clone.4', 'Bulk clone 4', outMalig$clone, fixed = T)
pdf(paste0(WD, 'endo_malig.pdf'), height=9, width=18)
ggplot(outMalig[c(seq(1,4), seq(9,16)),], aes(x = section, y = percentage, group = clone, shape = clone))+
    geom_point(size=8) +
    geom_line(size=3) +
    labs(x='Section', y='Percentage', title='Estimation of clone 4 cellular abundance', shape = 'Source of estimate')+
    scale_y_continuous(limit=c(0,60), breaks = seq(0, 60, 20)) +
    theme_classic() +
    oldham_theme() +
    theme(legend.spacing.y = unit(.5, 'cm')) +
    guides(shape = guide_legend(byrow = TRUE))
dev.off()
# }}}
# nonmalignant bulk
# {{{
me = data.frame(fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/Module_eigengenes_02-55-42.csv'))
bulkVecs = c('orange')
me = me[, c(1, which(colnames(me) %in% bulkVecs))]
me$Sample = as.numeric(gsub('.*_' , '', me$Sample))
snSlices = c(17, 53, 93, 117)
names(snSlices) = seq_along(snSlices)
inds = lapply(snSlices, function(slice) {
    order(abs(slice - me$Sample), decreasing=F)[1:4]
})
out = do.call(rbind, lapply(seq_along(inds), function(i) {
    x=inds[[i]]
    apply(me[x, ], 2, mean)
}))
out = data.frame(out)
nonmaligBulk = out
# }}}
# panel H, plot nonmalignant
# {{{
outNon = rbind(outP[seq(9,12), ], data.frame(source = rep('Bulk', 4), clone = rep('Endothelial', 4), percentage = nonmaligBulk$orange, section = c(17,53,93,117)))
cor(outNon$percentage[1:4], outNon$percentage[5:8])
outNon$clone = rep('Enodthelial, R=0.89',8)
outNon$percentage[1:4] = scale(outNon$percentage[1:4])
outNon$percentage[5:8] = scale(outNon$percentage[5:8])
outNon$section = factor(outNon$section, levels = c(17, 53, 93, 117))
outNon$source =  gsub('single_cell', 'Single-nucleus clone\n4:2:2, R=0.93', outNon$source)
outNon$source =  gsub('Bulk', 'Bulk endothelial', outNon$source)
pdf(paste0(WD, 'endo_nonmalig_legend.pdf'), height=9, width=18)
ggplot(outNon, aes(x =section, y=percentage, group = source, shape=source))+
    geom_point(size=8) +
    geom_line(size=3) +
    labs(x='Section', y='Scaled values', title='Estimation of endothelial cellular abundance', shape = 'Source of estimate')+
    scale_y_continuous(limit = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 1.5)) +
    theme_classic() +
    oldham_theme() +
    theme(legend.spacing.y = unit(.5, 'cm')) +
    guides(shape = guide_legend(byrow = TRUE))
dev.off()
# }}}
# }}}
# final plot will be showing the enrichments of the single nucleus clones, this time broken down fully
# {{{
# combine oligo
clust = gsub(11, 8, clust, fixed=T)
out = future_lapply(unique(clust), function(clu){
    ind = (clu == clust)
    future_apply(exprD, 1, function(x){
        wilcox.test(x[ind], x[!ind], alternative='greater')$p.value
    })
})
out = lapply(out, function(x) {
    names(x) = exprGene
    return(x)
})
out.bak =out
out = lapply(out, sort, decreasing = F)
out = lapply(out, function(x) {
    x = x[!(is.na(names(x)))]
    x[x<(0.01/(19477*11))]
})
source('~/code/git/GSEA_generic/enrichment_functions.R')
names(out) = seq(1, length(out))
outReal = lapply(out, names)
rm('mySets')
enrichs = enrichment_man(outReal, exprGene, '/home/shared/genesets/genesets_slim')
colnames(enrichs)[8:ncol(enrichs)] = names(sort(clustConv))
write.csv(enrichs, file = paste0(WD,'sn_enrich.csv'), row.names = F)

library('flashClust')
source("/home/patrick/code/git/GSEA_generic/GSEAfxsV3.r")
broadSets=getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
enrichs = enrichment_man_broad(outReal, exprGene, broadSets)
write.csv(enrichs, file = paste0(WD,'sn_enrich_broad.csv'), row.names = F)
# }}}
# final heatmap
# {{{
WD = '/home/patrick/code/git/SF10711_single_nucleus/'
enrichments = read.csv(paste0(WD, 'sn_enrich.csv'))
enrichments_broad = read.csv(paste0(WD, 'sn_enrich_broad.csv'))
colnames(enrichments_broad) = colnames(enrichments)
enrichments = rbind(enrichments, enrichments_broad)
keepSet = c('MOSET6799', 'MOSET6835', 'M39055', # clone 1
            'MOSET6892', 'M39052', 'MOSET6', # clone 2
            'MOSET7065', 'MOSET17', 'M1938', # clone 3
            'MOSET6935', 'MOSET7059', 'MOSET6962', # clone 4:1 and clone 4:2:1
            'MOSET6801', 'MOSET6858', 'MOSET6796', # clone 5
            'MOSET7', 'MOSET6942', 'MOSET7051', # astrocytes
            'MOSET7052', 'MOSET6936', 'MOSET8', # oligodendrocytes 1 and 2
            'MOSET6946', 'MOSET7053', 'MOSET6894', # microglia
            'MOSET9', 'MOSET7054', 'MOSET3', # neuron
            'MOSET6955', 'MOSET7058', 'MOSET6956' # endothelial
            )
nameSet = c('Verhaak: mesenchymal subtype', 'Tesar: OPC', 'Manno: radial glia',
            'Suva: astrocytoma program', 'Manno: OPC', 'Bachoo: astrocyte',
            'Monje: up in GBM line after NLGN3', 'Foster: mitochondrial proteins', 'Meissner: CpG promoters with\nhistone metyhylation',
            'Zeisel: ependymal', 'Kelley: ependymal', 'Gokce: ependymal',
            'Verhaak: proneural subtype', 'Noushmehr: up in proneural', 'Phillips: up in proneural',
            'Barres: astrocyte', 'Zeisel: astrocyte', 'Kelley: astrocyte',
            'Kelley: oligodendrocyte', 'Zeisel: oligodendrocyte', 'Barres: oligodendrocyte',
            'LaManno: microglia', 'Kelley: microglia', 'Suva: microglia',
            'Barres: neuron', 'Kelley: neuron', 'ABA: neuron',
            'HPA: endothelial', 'Kelley: endothelial', 'GTEX: endothelial'
            )
enrichments = enrichments[match(keepSet,enrichments$SetID),]
enrichments = as.matrix(enrichments[,8:ncol(enrichments)])
enrichments[enrichments == 0] = 1E-62
enrichments[enrichments < 1E-62] = 1E-62
enrichments = -log10(enrichments)
rownames(enrichments) = nameSet
colnames(enrichments) = gsub('([0-9])\\.([0-9])', '\\1:\\2', colnames(enrichments))
colnames(enrichments) = gsub('([a-z])\\.([0-9])', '\\1 \\2', colnames(enrichments))
colnames(enrichments) = gsub('(dendrocyte)\\ ([0-9])', '\\1:\\2', colnames(enrichments))
pdf(paste0(WD, 'single_nuc_enrichments.pdf'), height = 13, width = 10)
Heatmap(enrichments,
    name = '-log10 q-values',
    cluster_columns = FALSE,
    clustering_method_columns = 'ward.D',
    clustering_distance_columns = 'pearson',
    cluster_rows = FALSE,
    clustering_method_rows = 'ward.D',
    clustering_distance_rows = 'pearson',
    column_names_side = "top",
    column_names_centered = FALSE,
    column_names_gp = gpar(fontfamily = 'NimbusSan', cex = 2, col = c('#ff7f00', '#999999', '#dd95e8', '#4daf4a', '#90af41', '#a65729', '#e6ab02', '#e7298a', '#e729e6', '#1b9e77', '#701d02', '#1d4a4d')),
    row_title_gp = gpar(fontfamily = 'NimbusSan', cex = 2),
    row_names_gp = gpar(fontfamily = 'NimbusSan', cex = 1.5, col = rep(c('#ff7f00', '#999999', '#dd95e8', '#4daf4a', '#a65729', '#e6ab02', '#e7298a', '#1b9e77', '#701d02','#1d4a4d'), ea=3)),
    row_dend_width = unit(2, 'cm'),
    row_title = NULL,
    column_names_rot = 45, 
    row_split = rep(seq(1, 10), ea=3),
    rect_gp = gpar(col = 'black', lwd = 1),
    col = colorRamp2(breaks = c(0, 15, 66),
                            colors = c('#ffffff', '#fc9272', '#de2d26')
                    )
)
dev.off()
# }}}
# plot comparing kme to differential expression t-value
# {{{
# environment setup
# {{{
library('data.table')
library('ggplot2')
library('future')
library('future.apply')
library('qvalue')
library('dendextend')
library('RColorBrewer')
library('grid')
library('gridExtra')
library('limma')
library('ComplexHeatmap')
library('parallelDist')
library('kSamples')
library('ggsignif')
library('circlize')
library('phylogram')
library('slingshot')
library('igraph')
library('org.Hs.eg.db')
plan(multicore, workers=20)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
WD = '~/@patrick/SF10711/sn.rna.seq/figures/'
# }}}
# get FM clone genesets and nonmalignant genesets
# {{{
clones = c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
fmGenes = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
clones=c('red', 'violet', 'black', 'ivory', 'lightcyan')
fmGenes = lapply(clones, function(clone) fmGenes$Gene[which(fmGenes$'TopModPosFDR_0.12' == clone)])
names(fmGenes) = clones

fidGenes = lapply(list.files('/home/shared/genesets/fidelity_cell_types/hifi200', full.names=T), function(x) as.character(fread(x)[[1]]))
names(fidGenes) = gsub('_.*', '', list.files('/home/shared/genesets/fidelity_cell_types/hifi200'))
fidGenes= fidGenes[c(1, 10, 7, 9, 5)]
# we can instead get the genes from the dataset
clones = c('Astrocyte', 'Oligodendrocyte', 'Microglia', 'Neuron', 'Endothelial cell')
bulkGenes = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
clones=c('purple', 'tan', 'yellow', 'green', 'orange')
bulkGenes = lapply(clones, function(clone) bulkGenes$Gene[which(bulkGenes$'TopModPosBC_3.63e-08' == clone)])
names(bulkGenes) = clones
set1 = c(fmGenes, bulkGenes)
# }}}
# get expression matrix and cluster assignments
# {{{
dist = as.numeric(unlist(fread('~/@patrick/SF10711/sn.rna.seq/sanity/cell_cell_distance_with_errorbar_s2n_gt_1.txt')))
# read in expression matrix to determine wanted dimensions
expr = data.frame(fread('~/@patrick/SF10711/sn.rna.seq/sanity/log_transcription_quotients.txt'))
expr1 =expr[,1]
expr = expr[,-1]
expr=data.table(expr1, apply(expr, 2, function(x) exp(x))*1E6)
# format into distance matrix
zeros = seq(1,ncol(expr)-1)
delta = seq((ncol(expr)-3), 1)
indeces = c(1)
for(i in seq(1, ncol(expr)-3)){
    if(i==1){
        indeces [i+1] = indeces[i] + delta[i]
    } else{
        indeces [i+1] = indeces[i] + delta[i] +1
    }
}
indeces[1]=0
indeces[length(indeces)+1]=indeces[length(indeces)]+1
distMat=lapply(seq_along(indeces[-length(indeces)]), function(i) {
    out=dist[(indeces[i]+1) : indeces[i+1]]
    if(zeros[i] != 809){
        out=unlist(c(rep(0, zeros[i]),out))
    } else{
        out=unlist(c(rep(0, zeros[i]-1), out[2]))
    }
    return(out)
})
distMat=data.frame(do.call(cbind, distMat))
distMat=cbind(distMat, rep(0, nrow(distMat)))
# now need to mirror distance matrix
for(i in seq_len(nrow(distMat))) {
    for(j in seq_len(ncol(distMat))) {
        distMat[i, j] = distMat[j, i]
    }
}
distDf = as.dist(distMat)
# }}}
# harmonize gene names
# {{{
gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf', skip=5)
colnames(gtf)=c("chr", "db", "type", "start", "end", "score", "strand", "frame", "extra")
gtf=gtf[which(gtf$type=="gene"),]
gtf.ENSG=gsub(".*gene_id\\s\"\\\"*|\\..*", "", gtf$extra)
gtf.type=gsub(".*gene_type\\s\"\\\"*|\".*", "", gtf$extra)
gtf.name=gsub(".*gene_name\\s\"\\\"*|\".*", "", gtf$extra)
gtfOut=data.frame(gtf[,1], gtf.ENSG, gtf.type, gtf.name)
gtfOut$gtf.name= alias2SymbolTable(gtfOut$gtf.name, species = 'Hs')
gtfOut$gtf.name[is.na(gtfOut$gtf.name)] = gtfOut$gtf.ENSG[is.na(gtfOut$gtf.name)]
fmGenes=lapply(fmGenes, alias2SymbolTable, species='Hs')

expr1Gene = expr$expr1
expr = data.frame(gtfOut$gtf.name[match(as.character(unlist(expr$expr1)), gtfOut$gtf.ENSG)], expr)
colnames(expr)[1:2] = c('Gene', 'ENSG')
expr$Gene[is.na(expr$Gene)] = expr$ENSG[is.na(expr$Gene)]
exprC = fread('~/@patrick/SF10711/sn.rna.seq/190809_SF013_oldham_july_1/08_expression_matrix/expr.ensg.counts.csv')
exprC = exprC[match(expr$ENSG, exprC$Gene),]
exprC = data.frame(expr$Gene, exprC)
colnames(exprC)[1:2] = c('Gene', 'ENSG')

countThresh = future_apply(exprC[,-c(1,2)], 1, function(x) length(which(x<0.5)) / (ncol(exprC)-2))
threshInd = which(countThresh<0.9)
exprCD = exprC[threshInd,-c(1,2)]
exprD = expr[threshInd,-c(1,2)]
exprSum = future_apply(exprCD, 1, sum)
exprGene = expr$Gene[threshInd]
fmGenes = lapply(fmGenes, function(x) return(x[!(is.na(x))]))
# }}}
# get cluster definitions from hierarchical clustering
# {{{
attr(distDf, 'Labels') =paste0('x', seq(1, length(attr(distDf, 'Labels'))))
distHclust = hclust(distDf, method = 'ward.D')
dendP = as.dendrogram(distHclust)
clust = cutree(distHclust, k=12)
names(clust) = colnames(expr)[-c(1,2)]
# }}}

# input data
# {{{
bulkGenes = data.frame(fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv'))
bulkGenes$Gene= alias2SymbolTable(bulkGenes$Gene, species = 'Hs')
interGene = intersect(exprGene,bulkGenes$Gene)
bulkGenes = bulkGenes[match(interGene, bulkGenes$Gene),]
exprDT = exprD[match(interGene, exprGene),]
# 1) need vector of t-values for malignant cell clones
clust = cutree(distHclust, k=12)
out = future_lapply(unique(clust), function(clu){
    ind = (clu == clust)
    future_apply(exprDT, 1, function(x){
        t.test(x[ind], x[!ind], alternative='greater')$statistic
    })
})
freq = read.csv('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv')
rna=fread('/mnt/bdata/@patrick/SF10711/rna.seq/expression_matrices/Normalized_read_counts_using_RUVg_ERCC_K10Factors.csv', drop=seq(2,6))
colnames(rna)[1]='Gene'
rna = rna[match(interGene, rna$Gene), ]
rna = rna[, c(1, which(colnames(rna) %in% freq[,1])), with = F]
# 2) need vector of kME values for malignant cell clones
plotCors = list(
    astro = data.frame(bulk = bulkGenes$kMEpurple, singleNuc = unlist(out[[9]])),
    oligo = data.frame(bulk = bulkGenes$kMEtan, singleNuc = unlist(out[[11]])),
    micro = data.frame(bulk = bulkGenes$kMEyellow, singleNuc = unlist(out[[4]])),
    neuro = data.frame(bulk = bulkGenes$kMEgreen, singleNuc = unlist(out[[12]])),
    endo = data.frame(bulk = bulkGenes$kMEorange, singleNuc = unlist(out[[6]])),
    clone1 = data.frame(bulk = bulkGenes$kMEred, singleNuc = unlist(out[[2]])),
    clone2 = data.frame(bulk = bulkGenes$kMEviolet, singleNuc = unlist(out[[1]])),
    clone3 = data.frame(bulk = bulkGenes$kMEblack, singleNuc = unlist(out[[7]])),
    clone4 = data.frame(bulk = bulkGenes$kMEviolet, singleNuc = unlist(out[[5]])),
    clone5 = data.frame(bulk = bulkGenes$kMElightcyan, singleNuc = unlist(out[[10]]))
)
# alt, try different with correlation to clonal abundance
plotCors = list(
    clone1 = data.frame(bulk = as.numeric(cor(freq$clone.1, t(rna[,-1]))), singleNuc = unlist(out[[2]])),
    clone2 = data.frame(bulk = as.numeric(cor(freq$clone.2, t(rna[,-1]))), singleNuc = unlist(out[[1]])),
    clone3 = data.frame(bulk = as.numeric(cor(freq$clone.3, t(rna[,-1]))), singleNuc = unlist(out[[7]])),
    clone4 = data.frame(bulk = as.numeric(cor(freq$clone.4, t(rna[,-1]))), singleNuc = unlist(out[[5]])),
    clone5 = data.frame(bulk = as.numeric(cor(freq$clone.5, t(rna[,-1]))), singleNuc = unlist(out[[10]]))
)
xlabsCor = c('0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-95', '95+', '99+')
quants = c(seq(10, 90, 10), 95, 96,99)
# }}}
# {{{
# {{{
namesCor = c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Neurons', 'Endothelial cells', 'Clone 1', 'Clone 2', 'Clone 3', 'Clone 4', 'Clone 5')
colTitle = c('#ff7f00','#e6ab02','#999999','#e7298a','#dd95e8', '#1b9e77','#4daf4a','#701d02','#a65729','#1d4a4d')[c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9)]
gPlotcor = lapply(seq_along(plotCors), function(s) {
    t = plotCors[[s]]
    t$rankSingleNuc = rank(t$singleNuc)/nrow(t)*100
    t$quantRank = rep(NA, nrow(t))
    for(quant in quants){
        if(quant == 95){
            start = 90
            end = 95
        } else if(quant == 96){
            start = 95
            end = 100
        } else if(quant == 99){
            start = 99
            end =100
        } else{
            start = quant - 10
            end = quant
        }
        t$quantRank[t$rankSingleNuc >= start & t$rankSingleNuc <= end] = quant
    }
    t$quantRank = factor(t$quantRank, levels = quants)
    ggplot(t, aes(x = quantRank, y = bulk)) +
        geom_violin(scale='width', draw_quantiles = c(0.5)) + 
        labs(title = namesCor[s], x = '', y = '') +
        scale_y_continuous(breaks = seq(-1,1, 1), limits = c(-1,1))+
        scale_x_discrete(labels = xlabsCor)+
        theme_classic() +
        oldham_theme() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 30),
        plot.title = element_text(color = colTitle[s], size =60))
})
# }}}
# alt
# {{{
namesCor = c('Clone 1', 'Clone 2', 'Clone 3', 'Clone 4', 'Clone 5')
colTitle = c('#1b9e77','#4daf4a','#701d02','#a65729','#1d4a4d')[c(1, 3, 5, 2, 4)]
gPlotcor = lapply(seq_along(plotCors), function(s) {
    t = plotCors[[s]]
    t$rankSingleNuc = rank(t$singleNuc)/nrow(t)*100
    t$quantRank = rep(NA, nrow(t))
    for(quant in quants){
        if(quant == 95){
            start = 90
            end = 95
        } else if(quant == 96){
            start = 95
            end = 100
        } else if(quant == 99){
            start = 99
            end =100
        } else{
            start = quant - 10
            end = quant
        }
        t$quantRank[t$rankSingleNuc >= start & t$rankSingleNuc <= end] = quant
    }
    t$quantRank = factor(t$quantRank, levels = quants)
    ggplot(t, aes(x = quantRank, y = bulk)) +
        geom_violin(scale='width', draw_quantiles = c(0.5)) + 
        labs(title = namesCor[s], x = '', y = '') +
        scale_y_continuous(breaks = seq(-1,1, 1), limits = c(-1,1))+
        scale_x_discrete(labels = xlabsCor)+
        theme_classic() +
        oldham_theme() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 30),
        plot.title = element_text(color = colTitle[s], size =60))
})
# }}}
pdf(paste0(WD, 'violin_kme_cors_clon_abund_un.pdf'), height =13*1.3, width = 13*2.5)
plotL=list(grobs=gPlotcor, nrow=1, 
            bottom = textGrob('Single cell DE (t-value percentile)',  gp = gpar(fontsize = 65, fontface='bold', fontfamily='NimbusSan')),
            left = textGrob('Pseudobulk kME',  gp = gpar(fontsize = 65, fontface='bold', fontfamily='NimbusSan'), rot = 90)
)
margin = theme()
grid.arrange(do.call(arrangeGrob, plotL), padding = unit(c(0,0,0,0), "cm"))
dev.off()
# }}}
# }}}
