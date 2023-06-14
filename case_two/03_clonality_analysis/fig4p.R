setwd('~/@patrick/SF10711/figures')
library('stringr')
barOut=read.csv('tables/sf10711_cellular_abundance.csv')
colnames(barOut)=c('section', 'Nonmalignant', 'Clone 1', 'Clone 2', 'Clone 3', 'Clone 4', 'Clone 5')
barOutP=melt(barOut)
barOutP$section=factor(barOutP$section, levels=unique(barOutP$section))
oldham_theme=function() {theme(axis.line = element_line(colour = "black"),
	legend.title = element_blank(),
 	axis.text.x = element_text(size=40, color='black',  family='NimbusSan', margin=margin(t=10)),
 	axis.text.y = element_text(size=40, color='black', family='NimbusSan', margin=margin(r=10)),
 	axis.title.y = element_text(size=40, face='bold', family='NimbusSan', margin=margin(t=0, r=-10, b=0, l=0)),
 	axis.title.x = element_text(size=40, face='bold', family='NimbusSan', margin=margin(t=-10, r=0, b=0, l=0)),
 	plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan', margin=margin(t=-20, b=0)),
 	plot.subtitle = element_text(size=40,face="bold", hjust=.5 , family='NimbusSan', margin=margin(t=40, b=-40)),
 	axis.line.x = element_line(size=3),
 	axis.line.y = element_line(size=3),
 	plot.margin = unit(c(4, 2, 1, 2), "lines"),
 	legend.position="right",
	legend.key.size=unit(1.5, 'cm'),
 	legend.text=element_text(size=30, family='NimbusSan'))
}
barOutP$section=as.numeric(barOutP$section)
setwd('~/@patrick/SF10711/figures')
pdf('clonal_abundance.pdf', width=24, height=13)
print(ggplot(barOutP, aes(x=section, y=value, color=variable, fill=variable)) +
	theme_classic()+
	oldham_theme()+
	geom_bar(stat="identity", width=1, size=2)+
	scale_y_continuous(expand=c(0,0), breaks=c(0,1), labels=c(0,100))+
	scale_x_continuous(expand = c(0, 0), limits=c(1,85), breaks=c(1,85), labels=c(1,140))+
	labs(color='Clone',title = "Cellular fraction", x = "Section ID", y = "Percentage")+
	geom_vline(xintercept=46, size=2)+
	scale_fill_manual(values=c('white', '#ff7f00', '#984ea3', '#dd95e8', '#4daf4a', '#a65628'), name="Clones",labels=c('Nonmalignant', 'Clone 1', 'Clone 2', 'Clone 3', 'Clone 4', 'Clone 5'))+
	scale_color_manual(values=c('black', '#ff7f00', '#984ea3', '#dd95e8', '#4daf4a', '#a65628'), name="Clones",labels=c('Nonmalignant', 'Clone 1', 'Clone 2', 'Clone 3', 'Clone 4', 'Clone 5')))
dev.off()


setwd('~/bdata/SF10711/clonality/pyclone/outputs/tables')
loci=fread('cluster.tsv')
loci=reshape(loci[,-3], idvar='sample_id', timevar='cluster_id', direction='wide')
write.table(loci, row.names=F, quote=F, file='~/bdata/SF10711/figures/tables/pyclone_freq.tsv')

setwd('~/bdata/SF10711/clonality/pyclone_cITUP/outputs/tmp/tmp/tree/7')
tree=fread('results', skip='#clone_freq')    
colnames(tree)=c('c1', 'c2', 'c3', 'c4', 'c5') 
write.table(tree, row.names=F, quote=F, file='~/bdata/SF10711/figures/tables/citup_freq.tsv')
