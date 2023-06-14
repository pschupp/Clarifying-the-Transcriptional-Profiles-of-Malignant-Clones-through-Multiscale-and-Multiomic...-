library('data.table')
library('ggplot2')
WD = '~/@patrick/SF10711/figures/'
barOut=data.frame(fread('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv'))
# construct cumulative clonal abundances by addition
barOutIn=data.frame(barOut[,-c(1,2,3,4)], Malignant=1-barOut$nonmalignant, Clone.2=(barOut$clone.2+barOut$clone.3))
barOutIn=barOutIn[,c(4,5, 1,2,3)]
colnames(barOutIn)=c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')

corPlot=reshape2::melt(cor(barOutIn))
corPlot=data.table(corPlot, label=round(corPlot$value,2))
corPlot$Var1=gsub('\\.', ' ', corPlot$Var1)
corPlot$Var2=gsub('\\.', ' ', corPlot$Var2)
pdf(paste0(WD, 'cor_heatmap_cumulative.pdf'), height=13, width=13)
    print(ggplot(corPlot, aes(x=Var1, y=Var2, fill=value)) +
        theme_classic()+
        oldham_theme()+
        geom_tile() +
        geom_text(aes(Var1, Var2, label=label), color='black',size=7) +
        scale_fill_distiller(palette='RdBu', breaks=seq(-1,1,1), limits=c(-1,1)) +
        labs(title='Cumulative clonal\nabundance correlations', x='', y='', fill='Pearson\ncorrelation') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=40), 
              axis.text.y = element_text(size=40),
              plot.title = element_text(size=60),
#              legend.title=element_text(size=55, family='NimbusSan'), 
              legend.text=element_text(size=40, family='NimbusSan'))
)
dev.off()
