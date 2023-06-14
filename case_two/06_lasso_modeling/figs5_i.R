# read in enrichments 
# {{{
library('data.table')
library('ggplot2')
library('grid')
library('RColorBrewer')
library('gridExtra')
broadEnrich=fread('~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/stability_45/enrich_results_boot_model_fdr_threshold_45.csv')
ourEnrich=fread('~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/stability_45/enrich_results_broad_boot_model_fdr_threshold_45.csv')
cnvEnrich=fread('~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/stability_45/enrich_results_cnv_model_fdr_threshold_45.csv')
# }}}
# select genesets
# {{{
colnames(cnvEnrich)[1]='SetName'
cnvEnrich$SetName = gsub('\\.csv', '', cnvEnrich$SetName)
cnvEnrich$SetID=cnvEnrich$SetName
genesetL = list(
                Clone.1= c('Chr7 amplification','GOBP_TRANSLATIONAL_ELONGATION', 'Chr9 amplification','REACTOME_PTEN_REGULATION', 'Chr8 amplification', 'GOCC_MITOCHONDRION', 'ACEVEDO_LIVER_TUMOR_VS_NORMAL_ADJACENT_TISSUE_UP'),
                Clone.3 = c('GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_CHEMOTAXIS', 'GOBP_ENDOTHELIN_RECEPTOR_SIGNALING_PATHWAY','TURCAN_HYPOMETHYLATED_UPREGULATED_IN_IDH1MUT_CLINICAL_SAMPLES_236_GENES','Chr2q deletion', 'Chr18 deletion'),
                Clone.4 = c('kelley_Ependymal','DESCARTES_FETAL_LUNG_CILIATED_EPITHELIAL_CELLS',  'ZEISEL_EPENDYMAL'),
                Clone.5 = c('kelley_microglia', 'LaMANNO_MICROGLIA_MIDBRAIN', 'ZEISEL_MICROGLIA','Chr11q amplification'))
# subset to genesets of intrest
genesetL = list(
    Clone.1 = c('Chr7 amplification', 'Chr8 amplification', 'Chr9 amplification', 'MOSET7053','MOSET6799'), # ,'M39056','MOSET6946',
    Clone.2 = c('M15592', 'MOSET6899', 'MOSET6759', 'MOSET239'), 
    Clone.3 = c('Chr2q deletion', 'Chr18 deletion', 'MOSET7054','MOSET25','M535','M15709'), # 'MOSET9','MOSET17',
    Clone.4 = c('MOSET7059', 'MOSET6935', 'M40124'), #'M6695'
    Clone.5 = c('Chr11q amplification','MOSET6801','MOSET241','MOSET227', 'M10792')
)
ourEnrich = rbind(ourEnrich[which(ourEnrich$SetID %in% genesetL$Clone.1),-c(3,4,5,6,7)], ourEnrich[which(ourEnrich$SetID %in% unlist(genesetL[-1])),-c(3,4,5,6,7)])
broadEnrich = rbind(broadEnrich[which(broadEnrich$SetID %in% genesetL$Clone.1),-c(3,4,5,6,7)], broadEnrich[which(broadEnrich$SetID %in% unlist(genesetL[-1])),-c(3,4,5,6,7)])
enrich=rbind(ourEnrich, broadEnrich, cnvEnrich)
# function when going over these genesets to query whether positive or negative is higher and then take the max and record sign
enrichOut = lapply(seq_along(genesetL), function(i) {
    clone = names(genesetL)[i]
    sets = genesetL[[i]]
    temp = enrich[which(enrich$SetID %in% sets),]
    cloneInd = grep(toupper(clone), toupper(colnames(temp)))
    rows=temp[,1]
    temp = apply(temp[,-c(1,2)],2, as.numeric)
    minT = lapply(seq(1, ncol(temp),2), function(j){
        t1 = apply(temp[, j:(j+1)], 1, min)
        return(t1)
    })
    minT = do.call(cbind, minT)
    signT = lapply(seq(1, ncol(temp),2), function(j){
        t1 = apply(temp[, j:(j+1)], 1, which.min)
        return(t1)
    })
    signT = do.call(cbind, signT)
    signT[signT==1]='positive' 
    signT[signT==2]='negative' 
    minT=data.frame(rows, minT)
    colnames(minT)=c('SetID', unique(gsub('_.*', '', colnames(temp))))
    minT=reshape2::melt(minT)
    minT$value=-log10(minT$value)
    minT$value[reshape2::melt(signT)$value=='negative'] = minT$value[reshape2::melt(signT)$value=='negative'] * (-1)
    minT=minT[order(minT$value),]
    return(minT)
})
enrichOut = do.call(rbind, enrichOut)
genesetLUn=unlist(genesetL)
write.table(unique(enrichOut$SetID), file='~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/stability_45/lassoEnrichmentSetNames.csv', quote=F, row.names=F, col.names=F, sep=',')
enrichOut$SetName = factor(enrichOut$SetName, levels = genesetLUn)
# }}}
# plot enrichment plot
# {{{
enrichNameConv=read.csv('~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/stability_45/lassoEnrichmentSetNames.csv',header=F)
enrichNameConv$V2=gsub('\\\\n', '\n', enrichNameConv$V2)
enrichOut$SetName=enrichNameConv$V2[match(enrichOut$SetID, enrichNameConv$V1)]
enrichOrd = unlist(lapply(seq_along(genesetL), function(i) {
   temp=enrichOut[(enrichOut$variable == names(genesetL)[i]),]
   temp=temp[(temp$SetID %in% genesetL[[i]]),]
   return(unique(temp$SetName[order(temp$value, decreasing=T)]))
}))
enrichOrd=enrichOrd[seq(length(enrichOrd),1)]
enrichOut$SetName=factor(enrichOut$SetName, levels = enrichOrd)
# windsorize
enrichOut$value[enrichOut$value>30] = 30
enrichOut$value[enrichOut$value<(-30)] = (-30)
enrichOut$variable=gsub('\\.', ' ', enrichOut$variable)
temp = ggplot(enrichOut, aes(x = variable, y  = SetName, fill=value))+
          geom_tile() +
          theme_minimal()+
          oldham_theme()+
          scale_fill_distiller(palette = 'RdBu', limits=c(-30,30), breaks=seq(-30,30,30)) +
          guides(fill = guide_colourbar(title.position = "top",
                            title.hjust = .5,
                            label.position = "right"))+
          labs(x='', y='', title='Enrichments of LASSO model', fill='Q-value')+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=40), 
                legend.position = 'right',
                axis.text.y = element_text(size=40),
                plot.title = element_text(size=80, hjust = 1),
                legend.title=element_text(size=55, family='NimbusSan'), 
                legend.text=element_text(size=40, family='NimbusSan'),
                axis.title.x = element_text(size=55, face='bold', family='NimbusSan', margin=margin(t=0, r=0, b=0, l=0)),
                legend.key.size=unit(3, 'cm'))
grid.arrange(textGrob(''),
            textGrob ('Enrichments of LASSO model',
            gp = gpar(fontsize=60, fontfamily = 'NimbusSan',fontface="bold")),
            temp,
            heights = c(.07,0.0005, 1.1),
            padding=unit(5,'line'))
pdf('~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/stability_45/heatmap_lasso_model.pdf', height=20, width=20)
print(temp)
dev.off()
# }}}
