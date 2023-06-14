barOut=barOutIn=read.csv('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv')
barOutIn=data.frame(barOut[,-c(1,2,3,4)], Malignant=1-barOut$nonmalignant, Clone.2=(barOut$clone.2+barOut$clone.3))
barOutIn=barOutIn[,c(4,5,1,2,3)]
colnames(barOutIn)=c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
# Fisher's exact gene set enrichment
# {{{
# 1 read in enrichments
setwd('~/@patrick/SF10711/integration_analysis/network_deconvolution')
library('data.table')
library('forecast')
library('gridExtra')
library('grid')
library('RColorBrewer')
library('ggplot2')
library('egg')
library('readxl')

source('~/code/git/SF10711_clonal_integration_analysis/gseaplot2_oldham.R')
oldhamSets=list(fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/oldham_slim_fdr_02-55-42_FDR.csv'), 
    fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/oldham_slim_fdr_neg_02-55-42_FDR_NEG.csv'))
broadSets =list(data.table(read_excel('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/GSHyperG_BROAD_SETS_02-55-42_TOPMODPOSFDR.xlsx')), 
    fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/GSHyperG_BROAD_SETS_02-55-42_FDR_NEG.csv'))
kME=fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
# cnv sets
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

source('~/code/git/GSEA_generic/enrichment_functions.R')

cnvPos=lapply(cnvGenes, function(x) {
    lapply(c('red', 'violet', 'black', 'ivory', 'lightcyan'), function(mod){
        fisherTest(unique(x, na.rm=T), unique(kME$Gene[kME$'TopModPosFDR_0.12' %in% mod], na.rm=T), unique(kME$Gene, na.rm=T))
    })
})
cnvPos = do.call(cbind, cnvPos)
rownames(cnvPos) = c('red', 'violet', 'black', 'ivory', 'lightcyan')
cnvPos=data.frame(module = c('red', 'violet', 'black', 'ivory', 'lightcyan') ,cnvPos)
cnvNeg=lapply(cnvGenes, function(x) {
    lapply(c('red', 'violet', 'black', 'ivory', 'lightcyan'), function(mod){
        fisherTest(unique(x, na.rm=T), unique(kME$Gene[kME$'TopModNegFDR_0.12' %in% mod], na.rm=T), unique(kME$Gene, na.rm=T))
    })
})
cnvNeg = do.call(cbind, cnvNeg)
rownames(cnvNeg) = c('red', 'violet', 'black', 'ivory', 'lightcyan')
cnvNeg=data.frame(module = c('red', 'violet', 'black', 'ivory', 'lightcyan') ,cnvNeg)
# cnvNeg = melt(data.table(cnvNeg), id.var='module')

genesetL = list(red = c('Chr7.amplification', 'Chr8.amplification', 'Chr9.amplification','MOSET6946','MOSET7053','MOSET6799','M19062','M39056'),
               violet = c('MOSET239', 'MOSET6899', 'MOSET241', 'M18403', 'M25882', 'M10268'),
                black = c('Chr2q.deletion', 'Chr18.deletion', 'MOSET9','MOSET7054','MOSET6717','MOSET6879','M40158','M39049'),
                ivory = c('MOSET7059','MOSET6962','MOSET6961','MOSET6935','MOSET7041','M40071'),
                lightcyan = c('Chr11q.amplification', 'MOSET227','MOSET6801','MOSET244','M15206','M22036','M3009'))

# subset to genesets of intrest
oldhamSets = lapply(oldhamSets, function(x) x[which(x$SetID %in% unlist(genesetL)),c(1,2, match(names(genesetL), colnames(x))), with=F])
oldhamSets = lapply(oldhamSets, function(x) x[which(x$SetID %in% unlist(genesetL)),c(1,2, match(names(genesetL), colnames(x))), with=F])
oldhamSets[[1]] = rbind(oldhamSets[[1]], data.frame(SetID=colnames(cnvPos)[-1], SetName = colnames(cnvPos)[-1],t(cnvPos)[-1,]))
oldhamSets[[2]] = rbind(oldhamSets[[2]], data.frame(SetID=colnames(cnvNeg)[-1], SetName = colnames(cnvNeg)[-1],t(cnvNeg)[-1,]))
broadSets = lapply(broadSets, function(x) x[which(x$SetID %in% unlist(genesetL)),c(1,2, match(names(genesetL), colnames(x))), with=F])

# function when going over these genesets to query whether positive or negative is higher and then take the max and record sign
setEval = function(enrichL,geneS) {
    enrichL = lapply(enrichL, function(x) x[which(x$SetID %in% geneS),])
    enrichL[[2]] = enrichL[[2]][match(enrichL[[1]]$SetID, enrichL[[2]]$SetID), ]
    enrichL = lapply(enrichL, function(x) melt(x[,-1], id.var='SetName'))
    out = data.frame(matrix(nrow=nrow(enrichL[[1]]), ncol=2))
    for(i in seq_along(rownames(enrichL[[1]])) ) {
        ind = which.min(c(enrichL[[1]]$value[i], enrichL[[2]]$value[i]))
        out[i,1] = enrichL[[ind]]$value[i]
        if(ind == 1){
            out[i,2] = 'positive'
        }
        if(ind == 2){
            out[i,2] = 'negative'
        }
    }
    out=data.frame(enrichL[[1]][,-c(3)], out)
    return(out)
}

# create matrix with values and matrix with pos or neg
oldhamEval = lapply(genesetL, function(set) setEval(oldhamSets, set))
broadEval = lapply(genesetL, function(set) setEval(broadSets, set))
oldhamEval[[4]] = data.frame(oldhamEval[[4]])
outEval=lapply(seq_along(oldhamEval), function(i){rbind(oldhamEval[[i]],broadEval[[i]])})
names(outEval) = names(genesetL)
outEval = lapply(outEval, function(x) {
    x$variable=factor(x$variable, levels=c('red', 'violet', 'black', 'ivory', 'lightcyan'))
    return(x)
})
# write.table(unique(unlist(lapply(outEval, function(x) x$SetName))), file='~/@patrick/SF10711/integration_analysis/network_deconvolution/geneset_names.txt', quote=F, row.names=F, col.names=F, sep=',')
genesetNames = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/geneset_names.txt', header=F)
colnames(genesetNames)=c('SetName', 'updatedSetName')
genesetNames$updatedSetName = gsub('\\\\n', '\n', genesetNames$updatedSetName)
outEval = lapply(seq_along(outEval), function(i) {
    x=outEval[[i]]
    colnames(x) = c('GeneSet', 'Module', 'P.value', 'Sign')
    x$P.value = as.numeric(x$P.value)
    x$P.value=-log10(x$P.value)
    x$P.value[x$Sign=='negative'] = x$P.value[x$Sign=='negative'] * -1
    # windsorize for clarity
    x$P.value[x$P.value>30]=30
    x$P.value[x$P.value<(-30)]=-30
    x$GeneSet = genesetNames$updatedSetName[match(x$GeneSet, genesetNames$SetName)]
    x$GeneSet = factor(x$GeneSet, levels = x$GeneSet[x$Module==names(outEval)[i]][order(x$P.value[x$Module==names(outEval)[i]])])
    return(x)
})
pdf('null.pdf')
enrichPlotFisher = lapply(outEval, function(x){
    temp=ggplot(x, aes(x=Module, y=GeneSet, fill=P.value))+
              geom_tile() +
              theme_minimal()+
              oldham_theme()+
              scale_fill_distiller(palette = 'RdBu', limits=c(-30,30), breaks=seq(-30,30,30)) +
              guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))+
              labs(x='', y='', title='', fill='q-value')+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=40), 
                    legend.position = 'bottom',
                    axis.text.y = element_text(size=40),
                    plot.title = element_text(size=60, hjust = 1.2),
                    legend.title=element_text(size=55, family='NimbusSan'), 
                    legend.text=element_text(size=40, family='NimbusSan'),
 	                axis.title.x = element_text(size=55, face='bold', family='NimbusSan', margin=margin(t=0, r=0, b=0, l=0)),
                    legend.key.size=unit(3, 'cm'))

    grid.arrange(textGrob(''),
                textGrob ('Enrichment of module genes',
                gp = gpar(fontsize=60, fontfamily = 'NimbusSan',fontface="bold")),
                temp,
                heights = c(.07,0.0005, 1.1),
                padding=unit(5,'line'))
})
dev.off()
names(enrichPlotFisher) = c('red', 'violet', 'black', 'ivory', 'lightcyan')
# }}}
# load RNA
# {{{
rna=fread('/mnt/bdata/@patrick/SF10711/rna.seq/expression_matrices/Normalized_read_counts_using_RUVg_ERCC_K10Factors.csv', drop=seq(2,6))
colnames(rna)[1]='Gene'
meanExpr=apply(rna[,-seq(1,6)], 1, var)
rna=rna[order(meanExpr, decreasing=T),]
rnarSums=apply(rna[,-1], 1, sum)
rna=rna[-which(rnarSums<quantile(rnarSums,0.15)),]
rna=rna[-grep('^ERCC-', rna$Gene),]
rnarSumsN=apply(rna[,-1], 1, sum)
rna=as.data.frame(rna)
rna=rna[,c(1, which(colnames(rna) %in% barOut$ampseq.sample))]
rnaScale=t(future_apply(rna[,-c(1)], 1, log2))
# correct line below and make sure that the sections line up as expected
colnames(rnaScale)=gsub('SF10711)', '', colnames(rnaScale))
rnaScale=data.frame(Gene=rna$Gene, rnaScale)
# }}}
# plotSum function
# {{{
rnaPlot=data.frame(Gene=rnaScale$Gene, t(apply(rnaScale[,-1], 1, scale)))
plotSum = function(cloneC, meC, titleText){
    plot1=data.frame(xvar=seq(1,nrow(me)), ModuleEigengene=scale(me[,which(colnames(me)==meC), with=F]), ClonalAbundance=scale(cloneC))
    colnames(plot1)=c('xvar', 'ModuleEigengene', 'ClonalAbundance')
    plot2=reshape2::melt(plot1, id.var='xvar')
    plot2CloneMA=data.frame(xvar=plot2$xvar[plot2$variable=='ClonalAbundance'], 
                            value=ma(plot2$value[plot2$variable=='ClonalAbundance'], order=8), 
                            variable=rep('', nrow(plot1)))
    plot2MEMA=data.frame(xvar=plot2$xvar[plot2$variable=='ModuleEigengene'], 
                         value=ma(plot2$value[plot2$variable=='ModuleEigengene'], order=8), 
                         variable=rep('', nrow(plot1)))

    sub=paste0('Pearson correlation: ', 
                round(cor.test(plot1$ModuleEigengene, plot1$ClonalAbundance)$estimate,2), ' (p = ', 
                signif(cor.test(plot1$ModuleEigengene, plot1$ClonalAbundance)$p.value,2),')')
    if(meC=='ivory'){
        p1meC='#8B8B83' # darker colors from https://r-charts.com/colors/
    } else if(meC=='lightcyan'){
        p1meC='#00EEEE'
    } else{
        p1meC=meC
    }
    p1=ggplot(plot2, aes(x=xvar, y=value, color=variable))+
        geom_point(size=5) +
        theme_classic()+
        oldham_theme()+
        scale_x_continuous(breaks=c(1,nrow(plot1)), labels=c(1,140))+
        geom_line(data=plot2CloneMA, size=3, color='#6699CC')+
        geom_line(data=plot2MEMA, size=3, color=p1meC)+
        scale_color_manual(values=c(p1meC, '#6699CC'), labels=c('Module eigengene','Clonal abundance'))+
        labs(x='Section ID', 
             y='Z-scored values', 
             color='', 
             title='ME vs. clonal abundance', 
             subtitle=sub) +
        theme(legend.position = "bottom",
                axis.text.y = element_text(size=40),
                axis.text.x = element_text(size=40),
                axis.title.y = element_text(size=40),
                axis.title.x = element_text(size=40),
                plot.title = element_text(size=60),
                plot.subtitle = element_text(size=40),
                legend.text=element_text(size=40, family='NimbusSan'))


    kmeCol=paste0('kME', meC)

    plot3=t(rnaPlot[(rnaPlot$Gene %in% (kME$Gene[order(kME[,which(colnames(kME)==kmeCol), with=F], decreasing=T)[1:12]])),-1])
    colnames(plot3)=rnaScale$Gene[(rnaScale$Gene %in% (kME$Gene[order(kME[,which(colnames(kME)==kmeCol), with=F], decreasing=T)[1:12]]))]
    rownames(plot3)=seq(1,nrow(plot3))
    plot3=reshape2::melt(plot3)
    p3=ggplot(plot3, aes(x=Var1, y=value, color=Var2))+
        geom_line(size=2)+
        theme_classic()+
        oldham_theme()+
        scale_x_continuous(breaks=c(1,nrow(plot1)), labels=c(1, 140))+
        scale_color_manual(values=colorRampPalette(brewer.pal(9, 'Set1'))(12))+
        guides(color=guide_legend(ncol=4, byrow=TRUE))+
        labs(x='Section ID', y='Scaed expression', color='', title='Genes most highly\ncorrelated to ME')+
        theme(legend.position = "bottom",
                axis.text.y = element_text(size=40),
                axis.text.x = element_text(size=40),
                axis.title.y = element_text(size=40),
                axis.title.x = element_text(size=40),
                plot.title = element_text(size=60),
                legend.text=element_text(size=32, family='NimbusSan'))

    p4 = enrichPlotFisher[which(names(enrichPlotFisher)==meC)]
    plotName=paste0(meC, '_plot.pdf')
    pdf(plotName, height=25*.85, width=28.125*.85)
    grid.arrange(arrangeGrob(p1, p3, nrow=2),p4[[1]], ncol=2,widths=4:5, top=textGrob(as.character(titleText), gp=gpar(fontsize=80,fontface="bold", fontfamily='NimbusSan')))
    dev.off()
}
# }}}
# load module eigengenes
# {{{
me=fread('/mnt/bdata/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/Module_eigengenes_02-55-42.csv')
me2=fread('~/@patrick/SF10711/figures/fig6/old/alt/Bicor-None_signum0.33_minSize20_merge_ME_0.8_20246/Module_eigengenes_03-30-05.csv')
me$lightcyan=me2$lightcyan
# trim module eigenes to sections for which we clonal abundance data (amp-seq data)
me=me[which(me$Sample %in% barOut$ampseq.sample),]
# }}}
# plot lineplots
plotSum(barOutIn$Clone.1, 'red', 'Clone 1 is correlated with red ME')
plotSum(barOutIn$Clone.2, 'violet', 'Clone 2 is correlated with violet ME')
plotSum(barOutIn$Clone.3, 'black', 'Clone 3 is correlated with black ME')
plotSum(barOutIn$Clone.4, 'ivory', 'Clone 4 is correlated with ivory ME')
plotSum(barOutIn$Clone.5, 'lightcyan', 'Clone 5 is correlated with lightcyan ME')
