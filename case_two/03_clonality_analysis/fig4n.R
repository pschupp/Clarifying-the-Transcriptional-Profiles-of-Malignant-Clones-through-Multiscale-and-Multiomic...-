WD = '~/@patrick/SF10711/figures/'
library('ggplot2')
library('data.table')
setwd('~/@patrick/SF10711/cnv.analysis/rna.seq/cnvkit-rna/outputs')
i=1
for(fileN in list.files()[grep('*cns', list.files())]){
	file=fread(fileN)
	if(i==1){
		out=file[,c(1,2,3)]
	}
	out=cbind(out, file[,5])
	colnames(out)[(3+i)]=gsub('\\..*', '', fileN)
	i=i+1
}
sec=order(as.numeric(gsub('.*_', '', colnames(out)[grep('SF10711', colnames(out))])))+3
out=data.frame(out)[,c(1,2,3,sec)]
out$diff=out[,3]-out[,2]
out$change=paste('chr',out[,1],':', out[,2], '-', out[,3], '\nlength:', out$diff, sep='')
outPlot=reshape2::melt(out[,-c(1,2,3,100)])
cnsEx=data.frame(fread("SF10711_9_1_1.cns"))
# exome
facets=fread('~/@patrick/SF10711/figures/tables/facets_final.csv')
chr = 2
tempPlot=outPlot[which(gsub(':.*', '', gsub('chr', '', outPlot$change))==chr),]
tempCns=cnsEx[which(cnsEx$chromosome==chr),]
tempPlot = tempPlot[tempPlot$change == 'chr2:94811046-241817413\nlength:147006367',]
tempPlot$variable = as.numeric(gsub('.*_', '', tempPlot$variable))
facetsP = data.frame(variable = c(22,45,85,123), value = c(0.2482502, 0.1479187, 0.1, 0.100000))
pdf(paste0(WD, 'chr2qdel.pdf'), height=13*.4, width = 13*.7)
print(ggplot()  +
    geom_bar(data = facetsP, aes(variable, value), stat = 'identity', fill = 'grey', color='black', width = 5 ) +
    geom_point(data = tempPlot, aes(variable, value), size =3)+
    geom_smooth(data = tempPlot, aes(variable, value), size =3, color = 'black')+
    scale_y_continuous(expand = c(0,0), breaks = c(0, 0.3))+
    scale_x_continuous(expand = c(0,0), breaks=c(0, 140))+
    labs(x='Section ID', y='logR', title = '') + # title='Chr2q deletion estimation by exome and RNA-seq') +
    theme_classic() +
    oldham_theme()
    )
dev.off()
