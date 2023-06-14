library('data.table')
# load clonal abundance vectors
# {{{
barOut=barOutIn=read.csv('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv')
barOutIn=data.frame(barOut[,-c(1,2,3,4)], Malignant=1-barOut$nonmalignant, Clone.2=(barOut$clone.2+barOut$clone.3))
barOutIn=barOutIn[,c(4,5,1,2,3)]
colnames(barOutIn)=c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
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
library('future')
library('future.apply')
library('readxl')
source('~/code/git/SF10711_clonal_integration_analysis/gseaplot2_oldham.R')
oldhamSets=list(fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/oldham_slim_fdr_02-55-42_FDR.csv'), 
    fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/oldham_slim_fdr_neg_02-55-42_FDR_NEG.csv'))
broadSets =list(data.table(read_excel('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/GSHyperG_BROAD_SETS_02-55-42_TOPMODPOSFDR.xlsx')), 
    fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/GSHyperG_BROAD_SETS_02-55-42_FDR_NEG.csv'))
kME=fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
source('~/code/git/GSEA_generic/enrichment_functions.R')
# }}}
# subset to genesets of intrest
# {{{
genesetL = list(
                purple = c('M39074', 'MOSET7051', 'MOSET7', 'MOSET6942'),
                tan = c('MOSET7052', 'MOSET6936', 'M1714', 'M39037'),
                yellow = c('MOSET7053', 'MOSET6937', 'M39077', 'M39051'),
                green = c('MOSET9', 'MOSET7054', 'M2808', 'M735'),
                orange = c('MOSET7058', 'MOSET6955', 'MOSET6944', 'M39018')
                )
oldhamSets = lapply(oldhamSets, function(x) x[which(x$SetID %in% unlist(genesetL)),c(1,2, which(colnames(x) %in% unlist(c('red', 'violet', 'black', 'ivory', 'lightcyan', names(genesetL))))), with=F])
broadSets = lapply(broadSets, function(x) x[which(x$SetID %in% unlist(genesetL)), c(1,2, which(colnames(x) %in% unlist(c('red', 'violet', 'black', 'ivory', 'lightcyan', names(genesetL))))), with=F])
# }}}
# create dataframe for plotting
# {{{
# function when going over these genesets to query whether positive or negative is higher and then take the max and record sign
setEval = function(enrichL,geneS) {
    enrichL = lapply(enrichL, function(x) x[which(x$SetID %in% geneS),])
    enrichL[[2]] = enrichL[[2]][match(enrichL[[1]]$SetID, enrichL[[2]]$SetID), ]
    enrichL = lapply(enrichL, function(x) reshape2::melt(x[,-1], id.var='SetName'))
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

outEval=lapply(seq_along(oldhamEval), function(i){rbind(oldhamEval[[i]],broadEval[[i]])})
names(outEval) = names(genesetL)
genesetNames = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/geneset_names_nonmalig.txt', header=F, sep = ',', data.table = FALSE)
colnames(genesetNames)=c('SetName', 'updatedSetName')
genesetNames$updatedSetName = gsub('\\\\n', '\n', genesetNames$updatedSetName)

outEval = lapply(seq_along(outEval), function(i) {
    x=outEval[[i]]
    colnames(x) = c('GeneSet', 'Module', 'P.value', 'Sign')
    x$P.value = as.numeric(x$P.value)
    x$P.value=-log10(x$P.value)
    x$P.value[x$Sign=='negative'] = x$P.value[x$Sign=='negative'] * -1
    # windsorize
    x$P.value[x$P.value>30]=30
    x$P.value[x$P.value<(-30)]=-30
    x$GeneSet = genesetNames$updatedSetName[match(x$GeneSet, genesetNames$SetName)]
    x$GeneSet = factor(x$GeneSet, levels = x$GeneSet[x$Module==names(outEval)[i]][order(x$P.value[x$Module==names(outEval)[i]])])
    x$Module = factor(x$Module, levels = c('purple', 'tan', 'yellow', 'green', 'orange', 'red', 'violet', 'black', 'ivory', 'lightcyan'))
    print(i)
    return(x)
})

mods = c('purple', 'tan', 'yellow', 'green', 'orange')
for(i in seq_along(mods)){
    modsL = unlist(c(mods[i],c('red', 'violet', 'black', 'ivory', 'lightcyan')))
    outEval[[i]] = outEval[[i]][which(outEval[[i]]$Module %in% modsL),]
}
rnaPlot=data.frame(Gene=rnaScale$Gene, t(apply(rnaScale[,-1], 1, scale)))

enrichPlotFisher = lapply(outEval, function(x){
    temp = ggplot(x, aes(x=Module, y=GeneSet, fill=P.value))+
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
                    legend.key.size=unit(3, 'cm'))+
                theme(legend.position = 'none', plot.margin = margin(b=-35))
})
# }}}
# create line plot
# {{{
kmeCols=paste0('kME', mods)
plotLine = lapply(kmeCols, function(kmeCol){
    plot3=t(rnaPlot[(rnaPlot$Gene %in% (kME$Gene[order(kME[,which(colnames(kME)==kmeCol), with=F], decreasing=T)[1:12]])),-1])
    colnames(plot3)=rnaScale$Gene[(rnaScale$Gene %in% (kME$Gene[order(kME[,which(colnames(kME)==kmeCol), with=F], decreasing=T)[1:12]]))]
    rownames(plot3)=seq(1,nrow(plot3))
    plot3=reshape2::melt(plot3)
    p3=ggplot(plot3, aes(x=Var1, y=value, color=Var2))+
        geom_line(size=2)+
        theme_classic()+
        oldham_theme()+
        scale_x_continuous(breaks=c(1,nrow(plot3)), labels=c(1,81))+
        scale_color_manual(values=colorRampPalette(brewer.pal(9, 'Set1'))(12))+
        guides(color=guide_legend(ncol=4, byrow=TRUE))+
        labs(x='Section ID', y='Scaled expression', color='', title='')+
        theme(legend.position = "bottom",
                axis.text.y = element_text(size=40),
                axis.text.x = element_text(size=40),
                axis.title.y = element_text(size=40),
                axis.title.x = element_text(size=40),
                plot.title = element_text(size=60),
                legend.text=element_text(size=32, family='NimbusSan'))
})

pdf('purple_mod.pdf', height=13, width=30)
grid.arrange(plotLine[[1]], enrichPlotFisher[[1]],ncol=2, top = textGrob('Purple module snapshot', gp = gpar(fontsize=80, fontfamily = 'NimbusSan',fontface="bold")))
dev.off()

pdf('tan_mod.pdf', height=13, width=30)
grid.arrange(plotLine[[2]], enrichPlotFisher[[2]],ncol=2, top = textGrob('Tan module snapshot', gp = gpar(fontsize=80, fontfamily = 'NimbusSan',fontface="bold")))
dev.off()

pdf('yellow_mod.pdf', height=13, width=30)
grid.arrange(plotLine[[3]], enrichPlotFisher[[3]],ncol=2, top = textGrob('Yellow module snapshot', gp = gpar(fontsize=80, fontfamily = 'NimbusSan',fontface="bold")))
dev.off()

pdf('green_mod.pdf', height=13, width=30)
grid.arrange(plotLine[[4]], enrichPlotFisher[[4]],ncol=2, top = textGrob('Green module snapshot', gp = gpar(fontsize=80, fontfamily = 'NimbusSan',fontface="bold")))
dev.off()

pdf('orange_mod.pdf', height=13, width=30)
grid.arrange(plotLine[[5]], enrichPlotFisher[[5]],ncol=2, top = textGrob('Orange module snapshot', gp = gpar(fontsize=80, fontfamily = 'NimbusSan',fontface="bold")))
dev.off()

enrichPlotFisher[[2]], enrichPlotFisher[[3]],enrichPlotFisher[[4]], ncol=2)

p4 = enrichPlotFisher[which(names(enrichPlotFisher)==meC)]
plotName=paste0(meC, '_plot.pdf')
pdf(plotName, height=25, width=28.125)
 grid.arrange(arrangeGrob(p1, p3, nrow=2),arrangeGrob(p4[[1]], p4[[2]]), ncol=2,widths=4:5, top=textGrob(as.character(titleText), gp=gpar(fontsize=65,fontface="bold", fontfamily='NimbusSan')))
grid.arrange(arrangeGrob(p1, p3, nrow=2),p4[[1]], ncol=2,widths=4:5, top=textGrob(as.character(titleText), gp=gpar(fontsize=65,fontface="bold", fontfamily='NimbusSan')))
dev.off()
# }}}
