# dependencies
# {{{
library('data.table')
library('WGCNA')
library('gplots')
library('ggplot2')
library('reshape2')
library('ComplexHeatmap')
library('circlize')
library('RColorBrewer')
library('scales')
library('metagMisc')
library('stringr')
library('gridExtra')
oldham_theme=function() {theme(axis.line = element_line(colour = "black"),
	legend.title = element_blank(),
 	axis.text.x = element_text(size=40, color='black',  family='NimbusSan', margin=margin(t=10)),
 	axis.text.y = element_text(size=40, color='black', family='NimbusSan', margin=margin(r=10)),
 	axis.title.y = element_text(size=40, face='bold', family='NimbusSan', margin=margin(t=0, r=-20, b=0, l=0)),
 	axis.title.x = element_text(size=40, face='bold', family='NimbusSan', margin=margin(t=-10, r=0, b=0, l=0)),
 	plot.title = element_text(size=55,face="bold", hjust=.5, family='NimbusSan', margin=margin(t=0, b=0)),
 	plot.subtitle = element_text(size=40,face="bold", hjust=.5 , family='NimbusSan', margin=margin(t=0, b=0)),
 	axis.line.x = element_line(size=3),
 	axis.line.y = element_line(size=3),
 	plot.margin = unit(c(4, 2, 1, 2), "lines"),
 	legend.position="right",
	legend.key.size=unit(1.3, 'cm'),
 	legend.text=element_text(size=30, family='NimbusSan'))
}
# }}}
# read in exome mutation data
# {{{
setwd('~/@patrick/SF10711/figures')
Muts=read.csv("~/@patrick/SF10711/figures/tables/High_confidence_mutations_>200reads>10xfreq_in_blood.csv")
freqs=grep("_freq",colnames(Muts))
freqs=Muts[,freqs]
## Remove var.freq from the colnames:
colnames(freqs)=gsub("_freq","",colnames(freqs))
colnames(freqs)=gsub("_"," ",colnames(freqs),fixed=TRUE)
rownames(freqs)=Muts$sample
freqs=freqs[grep('SF10711_9', rownames(freqs)),]
freqs.cols.orig=colnames(freqs)
names=read.csv('~/@patrick/SF10711/figures/tables/mutation_anno.table.csv', sep='\t')
names.df=data.frame(names$Gene, names$protein)
names.p=unlist(lapply(as.character(names.df$names.protein), function(x) strsplit(x, ",")[[1]][1]))
names.df$names.protein=names.p
names.df$names.protein[which(is.na(names.df$names.protein))]="intronic"
names.df[78,2]="G245V" # correcting to more "canonical" AA change
colnames(freqs)=gsub("\\s.*", "", colnames(freqs))
names.df=names.df[match(colnames(freqs), names.df$names.Gene),]
colnames(freqs)=paste(names.df[,1], names.df[,2], sep=" ")
colnames(freqs)[grep("^NA", colnames(freqs))]=freqs.cols.orig[grep("^NA", colnames(freqs))]
cor.all=cor(freqs)
write.csv(cor.all, file="correlation.matrix.all.57.mutations.csv")
# }}}
# select mutations for presentation
# {{{
freqs.trim=freqs
freqs.trim=freqs.trim[, -grep("OR52H1|TSPAN13|OR11L1", colnames(freqs.trim))]
geneSelec=c('PTPRS', 'C10orf54', 'CARD14', 'GALNT2', 'THSD7B', 'LAMA2', 'MLX', 'PLEKHA5', 'MLL4', 'KRT17', 'OR8K3', 'INA', 'RYR2', 'CDKN1A', 'ATRX', 'TP53', 'TRIB2', 'TTN', 'TMCO4', 'PPIG', 'SSH1', 'PTPRZ1', 'ZSCAN10', 'HOXD4', 'IDH1', 'LRP2', 'RUFY1')
freqs.trim=freqs.trim[, grep(paste(geneSelec, collapse='|'), colnames(freqs.trim))]
conv=read.csv('~/@patrick/SF10711/figures/liftover_GRCh37_to_38/sf10711_muts_liftover.csv', row.names=NULL)
colnames(freqs.trim)=conv$new[match(gsub('\\s.*', '', colnames(freqs.trim)), conv$old)]
nameOut=colnames(freqs.trim)
cluster2=hclust(dist(1-cor(freqs.trim)))
clustDist=dist(1-cor(freqs.trim))
corMat2=cor(freqs.trim)
for(i in seq(1,nrow(corMat2))){
	corMat2[i,i]=1
}
groups1 = cutree(cluster2, k=5)
groups1=as.factor(groups1)
names(groups1)=gsub('\\.', ' ', names(groups1))
cluster.numb=length(unique(groups1))
Clusters=data.frame(paste("c", groups1, sep=""))
rownames(Clusters)=names(groups1)
colnames(Clusters)="Clusters"
map=list(Clusters=setNames(brewer.pal(cluster.numb, "Set1") , as.character(unique(Clusters$Clusters))))
colnames(corMat2)=rownames(corMat2)=names(groups1)
corMat2=corMat2[match(rownames(Clusters), rownames(corMat2)), match(rownames(Clusters), colnames(corMat2))]
breaksList = seq(-1, 1,by = 0.01)
colnames(corMat2)=gsub('\\.', ' ', colnames(corMat2))
# }}}
# create correlation heatmap
# {{{
colorText=rep('black', 23)
library('dendextend')
colorVec=brewer_pal(palette="Set1")(9)[c(4,7,5,6,3)]
colorVec[4]='#c1c400'
colFun=colorRamp2(seq(-1,1,length=7),rev(brewer.pal(n = 7, name = "RdBu")))
row_dend = color_branches(cluster2, k=5, col=colorVec)
pdf("high.vaf.dist.pearson.pdf", width=8, height=8)
ht=Heatmap(corMat2, col=colFun, rect_gp = gpar(lwd=3, col='black'), show_row_dend=F, column_labels=rep('', ncol(corMat2)), column_dend_height=unit(.75, 'cm'),
		cluster_columns = row_dend, column_split=5,row_order=cluster2$order, column_gap=unit(0,'mm'), column_dend_gp=gpar(lwd=5),
		clustering_distance_columns=clustDist, show_row_names=F,column_title = NULL,
		heatmap_legend_param = list(col_fun=colFun, title='Pearson correlation', direction = "horizontal", at=seq(-1,1,1),border=F, legend_width = unit(8, "cm"), title_gp=gpar(fontsize=20, fontface='bold', fontfamily='NimbusSan'), title_position = "topcenter", labels_gp=gpar(fontsize=20, fontfamily='NimbusSan')), 
		cell_fun = function(j, i, x, y, width, height, fill) {
			if((i==5)&(j==5))grid.text("1", x, y, gp = gpar(fontsize=40,fontfamily='NimbusSan', fontface='bold', col='white'))
			if((i==9)&(j==9))grid.text("2", x, y, gp = gpar(fontsize=40,fontfamily='NimbusSan', fontface='bold', col='white'))
			if((i==1)&(j==1))grid.text("3", x, y, gp = gpar(fontsize=40,fontfamily='NimbusSan', fontface='bold', col='white'))
			if((i==6)&(j==6))grid.text("4", x, y, gp = gpar(fontsize=40,fontfamily='NimbusSan', fontface='bold', col='white'))
			if((i==16)&(j==16))grid.text("5", x, y, gp = gpar(fontsize=40,fontfamily='NimbusSan', fontface='bold', col='white'))
		},
		top_annotation = HeatmapAnnotation(ta=anno_block(gp=gpar(fill=colorVec, lwd=5), labels=c('Cl. 1', 'Cl. 2', 'Cl.3', 'Cl. 4', 'Cl. 5'), labels_gp=gpar(fontsize=20, fontface='bold',fontfamily='NimbusSan'))))+
	rowAnnotation(ra=anno_text(rownames(corMat2), just="left",gp=gpar(col=colorText, fontfamily='NimbusSan', fontsize=14)))
draw(ht, heatmap_legend_side='bottom', column_title='Clustering of Amp-Seq VAFs', column_title_gp=gpar(fontsize=28, fontface = "bold", fontfamily='NimbusSan'))
dev.off()
# }}}
# create lineplots
# {{{
gnames=names(groups1)
groups1=paste0('c', groups1)
groups1=gsub('c2', 1, groups1)
groups1=gsub('c4', 2, groups1)
groups1=gsub('c1', 3, groups1)
groups1=gsub('c3', 4, groups1)
groups1=gsub('c5', 5, groups1)
groups1=as.numeric(groups1)
names(groups1)=gnames
titles=c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5')
orderVec=colnames(corMat2)[cluster2$order]
freqs.trim=freqs.trim[,c(match(orderVec, colnames(freqs.trim)))]
freqs.trim=data.frame(sample=as.numeric(gsub('.*_', '', rownames(freqs.trim))), freqs.trim)
colnames(freqs.trim)=gsub('\\.', ' ', colnames(freqs.trim))
lims=list(c(0,.4) ,c(0,.3),c(0,.9),c(0,.5) , c(0,.4))
spec_col=brewer_pal(palette='Set1')(9)[c(4,7,5,6,3)]
spec_col[4]='#c1c400'
subtitles = c(expression(paste('mean depth: 1.0x10'^'5 ', 'reads')), expression(paste('mean depth: 2.3x10'^'4 ', 'reads')), expression(paste('mean depth: 7.3x10'^'4 ', 'reads')), expression(paste('mean depth: 4.4x10'^'4 ', 'reads')), expression(paste('mean depth: 2.9x10'^'4 ', 'reads')))
pdf('lineplot_clusters.pdf', width=28, height=21)
i=1
p1=list()
p1[[1]]=1
p1[[2]]=2
p1[[3]]=3
p1[[4]]=4
p1[[5]]=5
for(ea in seq(1,5)){
	freqs.clust1=freqs.trim[,c(1,which(colnames(freqs.trim) %in% names(groups1[which(groups1==ea)])))]
	if(ncol(freqs.clust1)==2){
		freqs.clust1.m=data.frame(sample=freqs.clust1$sample, variable=rep(colnames(freqs.clust1)[2], nrow(freqs.clust1)), value=freqs.clust1[,2])
	} else{
		freqs.clust1.m=melt(freqs.clust1, id.vars='sample')
	}
	freqs.clust1.m[,1]=as.numeric(freqs.clust1.m[,1])
	freqs.clust1.m[,2]=as.factor(freqs.clust1.m[,2])
	freqs.clust1.m[,2]=factor(freqs.clust1.m[,2], levels=as.character(unique(freqs.clust1.m[,2])))
	p1[[i]]=(ggplot(freqs.clust1.m, aes(x=sample, y=value, group=variable,color=variable)) +
		theme_classic()+
		oldham_theme()+
		theme(plot.title=element_text(hjust=1))+
#		guides(colour=guide_legend(ncol=3, byrow=T))+
		geom_line(size=2)+
		labs(x="Section ID", y="VAF", title=titles[i], subtitle = subtitles[ea])+
		scale_x_continuous(expand = c(0, 0), limits=c(4,129), breaks=c(1,129), labels=c(1,140))+
		scale_y_continuous(expand = c(0, 0), limits=lims[[i]], breaks=lims[[i]], labels=lims[[i]])+
		theme(plot.title=element_text(color=spec_col[i]))+
		scale_color_brewer(palette="Set1"))
	i=i+1
}
# }}}
# create figure 2m
# {{{
freqs.2i=freqs.trim[,grep("ATRX|TP53|IDH1", colnames(freqs.trim))]
freqs.2i=data.frame(sections=seq(1,85), freqs.2i, ATRX.div.2=freqs.2i[,1]/2, TP53.div.2=freqs.2i[,2]/2, PredictedChr2qloss=(apply(data.frame(freqs.2i[,1], freqs.2i[,2]), 1, mean)/2) - freqs.2i[,3])
freq2i.plot=melt(freqs.2i[,which(colnames(freqs.2i) %in% c('sections', 'IDH1.R132H', 'ATRX.div.2', 'TP53.div.2', 'PredictedChr2qloss'))], id.vars=c(1), measure.vars=c(2:5))
freq2i.plot[,1]=as.numeric(freq2i.plot[,1])
freq2i.plot$variable=gsub("IDH1.R132H", "IDH1 R132H\n(Chr2) [cluster 3]", freq2i.plot$variable)
freq2i.plot$variable=gsub("ATRX.div.2", "ATRX intronic/2\n(ChrX) [cluster 1]", freq2i.plot$variable)
freq2i.plot$variable=gsub("TP53.div.2", "TP53 G245V/2\n(Chr17) [cluster 1]", freq2i.plot$variable)
freq2i.plot$variable=gsub("PredictedChr2qloss", "Discordance\n(CNV related)", freq2i.plot$variable)
freq2i.plot$variable=factor(freq2i.plot$variable, levels=c(unique(freq2i.plot$variable)[2], unique(freq2i.plot$variable)[3], unique(freq2i.plot$variable)[1], unique(freq2i.plot$variable)[4]))
write.table(freqs.2i, file='chr2qloss.tsv', sep='\t')
p1[[6]]=(ggplot(freq2i.plot, aes(x=sections, y=value, group=variable,color=variable)) +
	theme_classic()+
	oldham_theme()+
	theme(plot.title=element_text(margin=margin(t=-20, b=10)), legend.position='right')+
	guides(colour=guide_legend(ncol=1, byrow=T))+
	geom_line(size=3) + 
	labs(x="Section ID", y="VAF", title='Discordance of\ntrunkal mutations') +
	scale_x_continuous(expand = c(0, 0), limits=c(1,nrow(freqs.2i)), breaks=c(2,nrow(freqs.2i)), labels=c(1,140)) + 
	scale_y_continuous(expand = c(0, 0), limits=c(0,0.42), breaks=c(0,0.4))+
	geom_vline(xintercept=14)+
	geom_vline(xintercept=30)+
	geom_vline(xintercept=56)+
	geom_vline(xintercept=83)+
	scale_color_manual(values=c('#fb9a99', '#e31a1c', "#33a02c", "#a65628"))
)
grid.arrange(	p1[[1]]+theme(plot.margin = unit(c(1,1,0,0.5), "cm")), 
				p1[[2]]+theme(plot.margin = unit(c(1,1,0,0.5), "cm")), 
				p1[[3]]+theme(plot.margin = unit(c(1,1,0,0.5), "cm")), 
				p1[[4]]+theme(plot.margin = unit(c(1,1,0,0.5), "cm")),
				p1[[5]]+theme(plot.margin = unit(c(1,1,0,0.5), "cm")), 
				p1[[6]]+theme(plot.margin = unit(c(1,0,0,0.5), "cm")), ncol = 2, padding=unit(0, 'line'))
# marginMin=theme(plot.margin = unit(c(1,1,0,0.5), "cm"))
dev.off()
# }}}
