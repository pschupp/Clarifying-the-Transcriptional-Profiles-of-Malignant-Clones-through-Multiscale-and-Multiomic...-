library('data.table')
library('ggplot2')
library('stringr')
library('gridExtra')
library('RColorBrewer')
source('/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R')
WD = '~/@patrick/SF10711/figures/'
cnvkit=fread('~/@patrick/SF10711/figures/tables/cnvkit_final.csv')
cnvkit=cnvkit[-5, c(1,grep('avg', colnames(cnvkit))), with=F]
facets=fread('~/@patrick/SF10711/figures/tables/facets_final.csv')
facets=facets[-5,-c(2,3), with=F]
plot=data.frame(melt(cnvkit, id.vars='chromosome'), melt(facets, id.vars='chromosome')[,-c(1,2)])
colnames(plot)=c('Var2.1','slice' ,'value.1', 'value')
plot$Var2.1=as.factor(plot$Var2.1)
oldham_theme=function() {theme(axis.line = element_line(colour = "black"),
	legend.title = element_blank(),
 	axis.text.x = element_text(size=40, color='black',  family='NimbusSan', margin=margin(t=10)),
 	axis.text.y = element_text(size=40, color='black', family='NimbusSan', margin=margin(r=10)),
 	axis.title.y = element_text(size=40, face='bold', family='NimbusSan', margin=margin(t=0, r=-10, b=0, l=0)),
 	axis.title.x = element_text(size=40, face='bold', family='NimbusSan', margin=margin(t=-10, r=0, b=0, l=0)),
 	plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan', margin=margin(t=-20, b=-40)),
 	plot.subtitle = element_text(size=40,face="bold", hjust=.5 , family='NimbusSan', margin=margin(t=40, b=-40)),
 	axis.line.x = element_line(size=3),
 	axis.line.y = element_line(size=3),
 	plot.margin = unit(c(0, 2, 1, 2), "lines"),
 	legend.position="right",
	legend.key.size=unit(1.5, 'cm'),
 	legend.text=element_text(size=30, family='NimbusSan'))
}
pdf(paste0(WD, 'exome_meth.pdf'), width=13, height=13)
titleP=ggplot()+
	ggtitle("Concordance of exome")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan'))
titleP2=ggplot()+
	ggtitle("and RNA-seq CNV calls")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan'))
titleP3=ggplot()+
	ggtitle("N=4 sections, Pearson's r=0.97")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=40,face="bold", hjust=.5, family='NimbusSan', margin=margin(b=-40)))
mainP=ggplot(plot, aes(x=value.1, y=value, color=Var2.1))+
	geom_point(size=10)+
	theme_classic()+
	oldham_theme()+
	scale_color_manual(labels=c('Chr2q del', 'Chr7 amp', 'Chr8q Amp', 'Chr9 Amp', 'Chr11q Del', 'Chr18p Del'),values=brewer.pal(7,'Set1')[-6])+
	scale_y_continuous(expand=c(0,0),limits=c(-0.02,1),breaks=c(0,1))+
	scale_x_continuous(expand=c(0,0),limits=c(0,1),breaks=c(1))+
	labs(x = "Mean frequency (RNA-seq)", y = "Mean frequency (Exome)")+ 
    theme(
            legend.position="top",
            legend.title = element_blank(),
            legend.spacing.x = unit(.5, 'cm'),
            legend.key.size=unit(2, 'cm'),legend.box="vertical", legend.margin=margin()
        ) + 
        guides(color=guide_legend(nrow=3,byrow=TRUE, ncol =3))
grid.arrange(titleP, titleP2, titleP3, mainP,ncol = 1, heights = c(1.2/20, 1.2/20, 0.9/20, 16.7/20))
dev.off()
pdf(paste0(WD, 'exome_meth_no_leg.pdf'), width=13, height=13)
titleP=ggplot()+
	ggtitle("Concordance of exome")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan'))
titleP2=ggplot()+
	ggtitle("and RNA-seq CNV calls")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan'))
titleP3=ggplot()+
	ggtitle("N=4 sections, Pearson's r=0.97")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=40,face="bold", hjust=.5, family='NimbusSan', margin=margin(b=-40)))
mainP=ggplot(plot, aes(x=value.1, y=value, color=Var2.1))+
	geom_point(size=10)+
	theme_classic()+
	oldham_theme()+
	scale_color_manual(labels=c('Chr2q del', 'Chr7 amp', 'Chr8q Amp', 'Chr9 Amp', 'Chr11q Del', 'Chr18p Del'),values=brewer.pal(7,'Set1')[-6])+
	scale_y_continuous(expand=c(0,0),limits=c(-0.02,1),breaks=c(0,1))+
	scale_x_continuous(expand=c(0,0),limits=c(0,1),breaks=c(1))+
	labs(x = "Mean frequency (RNA-seq)", y = "Mean frequency (Exome)")+ 
    theme(
            legend.position="none",
            legend.title = element_blank(),
            legend.spacing.x = unit(.5, 'cm'),
            legend.key.size=unit(2, 'cm'),legend.box="vertical", legend.margin=margin()
        ) + 
        guides(color=guide_legend(nrow=3,byrow=TRUE, ncol =3))
grid.arrange(titleP, titleP2, titleP3, mainP,ncol = 1, heights = c(1.2/20, 1.2/20, 0.9/20, 16.7/20))
dev.off()


ea=7
setwd(paste0('~/@patrick/SF10711/clonality/pyclone_cITUP/outputs/tmp/tmp/tree/', ea))
tree=fread('results', skip='#clone_freq', fill = TRUE)
colnames(tree)=c('c0', 'c1', 'c2', 'c3', 'c4')
out=data.frame(freqs.2i$PredictedChr2qloss, tree$c2)
colnames(out)=c('ampseq', 'methylation')
cor(chr17ampVec, meth17)
pdf(paste0(WD, 'chr17_amp_meth.pdf') ,width=13,height=13)
titleP=ggplot()+
	ggtitle("Concordance of amp-seq and")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan'))
titleP2=ggplot()+
	ggtitle("RNA-seq chr2q deletion calls")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan', margin=margin(b=-40)))
titleP3=ggplot()+
	ggtitle("N=85 sections, Pearson's r=0.98")+
	geom_point() +
	theme_void() +
 	theme(plot.title = element_text(size=40,face="bold", hjust=.5, family='NimbusSan', margin=margin(b=-40)))
mainP=ggplot(out, aes(y=ampseq, x=methylation))+
	geom_point(size=10)+
	theme_classic()+
	oldham_theme()+
	scale_y_continuous(expand=c(0,0),limits=c(0,0.2),breaks=c(0,0.2))+
	scale_x_continuous(expand=c(0,0),limits=c(0,0.2),breaks=c(0.2))+
	labs(x = "Mean frequency (RNA-seq)", y = "Mean frequency (Amp-seq)")
grid.arrange(titleP, titleP2, titleP3, mainP,ncol = 1, heights = c(1.2/20, 1.2/20, 0.9/20, 16.7/20))
dev.off()


setwd('~/bdata/SF10711/figures')
load('/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/SF10711/downsample_coverage_freq_error/freqErrorCoverage_TP53_100.RData')
errPerCov=errPerCov[order(as.numeric(gsub("coverage_|x", "", names(errPerCov))))]
names(errPerCov)=gsub("coverage_", "", names(errPerCov))
df=melt(errPerCov)    
colnames(df) <- c("RMSE", "Coverage")                           
df$Coverage <- factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])     
df$RMSE=as.numeric(df$RMSE)
pdf('resampleTP53.pdf', width=17, height=13)
print(ggplot(df, aes(x=Coverage, y=RMSE, group=Coverage))+
	theme_classic()+
	oldham_theme()+
	geom_boxplot(fill ='white', color = "black", notch = TRUE, outlier.alpha=0, lwd=1.5)+
	scale_y_continuous(expand = c(0, 0), limits=c(0,0.08), breaks=c(0,0.08))+
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family='NimbusSan'), 
	plot.title = element_text(margin=margin(t=0, b=-40)),
 	plot.subtitle = element_text(margin=margin(t=60, b=10)),
)+
	labs(title='TP53 VAF RMSE between\nfull and downsampled coverage', subtitle = "1000 resamples per coverage"))
dev.off()
load('/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/SF10711/downsample_coverage_freq_error/freqErrorCoverage_IDH1_100_resamples.RData')
errPerCov=errPerCov[order(as.numeric(gsub("coverage_|x", "", names(errPerCov))))]
names(errPerCov)=gsub("coverage_", "", names(errPerCov))
df=melt(errPerCov)    
colnames(df) <- c("RMSE", "Coverage")                           
df$Coverage <- factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])     
df$RMSE=as.numeric(df$RMSE)
pdf('resampleIDH1.pdf', width=17, height=13)
print(ggplot(df, aes(x=Coverage, y=RMSE, group=Coverage))+
	theme_classic()+
	oldham_theme()+
	geom_boxplot(fill ='white', color = "black", notch = TRUE, outlier.alpha=0, lwd=1.5)+
	scale_y_continuous(expand = c(0, 0), limits=c(0,0.08), breaks=c(0,0.08))+
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family='NimbusSan'), 
	plot.title = element_text(margin=margin(t=0, b=-40)),
 	plot.subtitle = element_text(margin=margin(t=60, b=10)),
)+
	labs(title='IDH1 VAF RMSE between\nfull and downsampled coverage', subtitle = "1000 resamples per coverage"))
dev.off()
