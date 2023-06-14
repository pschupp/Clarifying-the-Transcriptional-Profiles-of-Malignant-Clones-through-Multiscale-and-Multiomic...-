# downsampling analysis
# {{{
WD = '~/@patrick/SF10711/figures/'
library('reshape2')
library('ggplot2')
library('gridExtra')
# TP53 RMSE plot
# {{{
base::load('/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF10711/downsample_coverage_freq_error/freqErrorCoverage_TP53_100.RData')
errPerCov=errPerCov[order(as.numeric(gsub("coverage_|x", "", names(errPerCov))))]
names(errPerCov)=gsub("coverage_", "", names(errPerCov))
df=melt(errPerCov)    
colnames(df) <- c("RMSE", "Coverage")                           
df = df[!(df[,2] %in% c('10000x', '11000x', '12000x')),]
df$Coverage <- factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])     
df$RMSE=as.numeric(df$RMSE)
tp53RMSE = ggplot(df, aes(x=Coverage, y=RMSE, group=Coverage))+
	theme_classic()+
	oldham_theme()+
	geom_boxplot(fill ='white', color = "black", notch = TRUE, outlier.alpha=0, lwd=1.5)+
	scale_y_continuous(expand = c(0, 0), limits=c(0,0.08), breaks=c(0,0.08))+
	theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            plot.margin = margin(r = 4, l = 25, t = 0, b = 15),
         	plot.title = element_text(margin=margin(t=0, b=-40)),
 	        plot.subtitle = element_text(margin=margin(t=60, b=10)))+
	labs(title='Effect on TP53 VAF between\nfull and downsampled coverage', subtitle = "1000 resamples per coverage")
# }}}
# TP53 correlation plot
# {{{
base::load('/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF10711/downsample_coverage_freq_correlation/freqCorSelfCov_TP53_1000_resamples.RData', verbose =TRUE)
corPerCoverage=corPerCoverage[order(as.numeric(gsub("coverage_|x", "", names(corPerCoverage))))]
names(corPerCoverage)=paste0(gsub("coverage_", "", names(corPerCoverage)), 'x')
df=melt(corPerCoverage)    
df = df[!(df[,2] %in% c('6000x')),]
colnames(df) = c("Correlation", "Coverage") 
df$Coverage = factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])     
df$Correlation = as.numeric(df$Correlation)
tp53Corr = ggplot(df, aes(x=Coverage, y=Correlation, group=Coverage))+
	theme_classic()+
	oldham_theme()+
	geom_boxplot(fill ='white', color = "black", notch = TRUE, outlier.alpha=0, lwd=1.5)+
	scale_y_continuous(expand = c(0, 0), limits=c(-0.2,1), breaks=c(0, 1))+
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family='NimbusSan'), 
            plot.margin = margin(l = 30, t = 40),
            axis.title.y = element_text(margin = margin(t = 0, r = 60, b = 0, l = 0)))
# }}}
pdf(paste0(WD, 'tp53_downsample.pdf'), width = 13, height = 13)
grid.arrange(tp53RMSE, tp53Corr, ncol = 1)
dev.off()
# IDH1 RMSE plot
# {{{
base::load('/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF10711/downsample_coverage_freq_error/freqErrorCoverage_IDH1_100_resamples.RData')
errPerCov=errPerCov[order(as.numeric(gsub("coverage_|x", "", names(errPerCov))))]
names(errPerCov)=gsub("coverage_", "", names(errPerCov))
df=melt(errPerCov)    
colnames(df) <- c("RMSE", "Coverage")
df = df[!(df[,2] %in% c('11000x', '12000x')),]
df$Coverage <- factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])     
df$RMSE=as.numeric(df$RMSE)
idh1RMSE = ggplot(df, aes(x=Coverage, y=RMSE, group=Coverage))+
	theme_classic()+
	oldham_theme()+
	geom_boxplot(fill ='white', color = "black", notch = TRUE, outlier.alpha=0, lwd=1.5)+
	scale_y_continuous(expand = c(0, 0), limits=c(0,0.08), breaks=c(0,0.08))+
	theme(axis.title.x = element_blank(), 
            plot.margin = margin(r = 4,l = 25, t = 0, b = 15),
            axis.text.x = element_blank(), 
         	plot.title = element_text(margin=margin(t=0, b=-40)),
 	        plot.subtitle = element_text(margin=margin(t=60, b=10)))+
	labs(title='Effect on IDH1 VAF between\nfull and downsampled coverage', subtitle = "1000 resamples per coverage")
# }}}
# IDH1 correlation plot
# {{{
base::load('/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF10711/downsample_coverage_freq_correlation/freqCorSelfCov_IDH1_1000_resamples.RData', verbose =TRUE)
corPerCoverage=corPerCoverage[order(as.numeric(gsub("coverage_|x", "", names(corPerCoverage))))]
names(corPerCoverage)=paste0(gsub("coverage_", "", names(corPerCoverage)), 'x')
df=melt(corPerCoverage)    
colnames(df) = c("Correlation", "Coverage")
df$Coverage = factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))]) 
df$Correlation = as.numeric(df$Correlation)
idh1Corr = ggplot(df, aes(x=Coverage, y=Correlation, group=Coverage))+
	theme_classic()+
	oldham_theme()+
	geom_boxplot(fill ='white', color = "black", notch = TRUE, outlier.alpha=0, lwd=1.5)+
	scale_y_continuous(expand = c(0, 0), limits=c(0,1), breaks=c(0,1))+
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family='NimbusSan'), 
            plot.margin = margin(l = 30, t = 40),
            axis.title.y = element_text(margin = margin(t = 0, r = 60, b = 0, l = 0)))
# }}}
pdf(paste0(WD, 'idh1_downsample.pdf'), width = 13, height = 13)
grid.arrange(idh1RMSE, idh1Corr, ncol = 1)
dev.off()
