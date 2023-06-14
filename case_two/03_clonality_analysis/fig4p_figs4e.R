# FACETS vectors of CNA
# {{{
# read in imputed CNVS from FACETS via exome data
cnv22=read.csv('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-1-22.cval.450.csv')
cnv46=read.csv('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-1-46.cval.450.csv')
cnv85=read.csv('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-2-85.cval.500.csv')
cnv123=read.csv('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-2-123.cval.450.csv')
# processing unique cnvs
cnv22$chrom=gsub(23, 'X', cnv22$chrom)
cnv46$chrom=gsub(23, 'X', cnv46$chrom)
cnv85$chrom=gsub(23, 'X', cnv85$chrom)
cnv123$chrom=gsub(23, 'X', cnv123$chrom)
# unique cnv parsing
temp1=cnv22[which(cnv22$Minor.Copy.Number!=1|cnv22$Major.Copy.Number!=1),]
cnvs=data.frame(temp1,sec=rep('s22', nrow(temp1)))
temp1=cnv46[which(cnv46$Minor.Copy.Number!=1|cnv46$Major.Copy.Number!=1),]
cnvs=rbind(cnvs,data.frame(temp1,sec=rep('s46', nrow(temp1))))
temp1=cnv85[which(cnv85$Minor.Copy.Number!=1|cnv85$Major.Copy.Number!=1),]
cnvs=rbind(cnvs,data.frame(temp1,sec=rep('s85', nrow(temp1))))
temp1=cnv123[which(cnv123$Minor.Copy.Number!=1|cnv123$Major.Copy.Number!=1),]
cnvs=rbind(cnvs,data.frame(temp1,sec=rep('s123', nrow(temp1))))
cnvs=cnvs[order(cnvs$chrom, cnvs$start),]
cnvs=rbind(cnvs, rep(0, ncol(cnvs)))
cnvKit=fread('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/outputs/called_cnvs.csv')
cnvKit=cnvKit[,-1]
cnvs=cnvs[which(cnvs$chrom %in% cnvKit$chromosome),]
cnvs$chrom=as.numeric(cnvs$chrom)
cnvs$sec=as.numeric(gsub('s', '', cnvs$sec))
cnvs=cnvs[,c(1,seq(10,16))]
cnvs_u=data.frame(matrix(nrow=1, ncol=ncol(cnvs)))
colnames(cnvs_u)=colnames(cnvs)
closeness=0.05# snps with start AND stop positions within X% will be merged with each other 
for(i in unique(cnvs$chrom)){
	inTemp=cnvs[which(cnvs$chrom==i),]
	inVec=which((inTemp$start>=((1-closeness)*inTemp$start[1]))&(inTemp$end<=((1+closeness)*inTemp$end[1])))
	if(length(inVec)==0){inVec=1}
	while(nrow(inTemp)!=0){
		temp=inTemp[inVec,]
		temp$start=round(mean(temp$start))
		temp$end=round(mean(temp$end))
		cnvs_u=rbind(cnvs_u,temp)
		inTemp=inTemp[-inVec,]
		inVec=which((inTemp$start>=((1-closeness)*inTemp$start[1]))&(inTemp$end<=((1+closeness)*inTemp$end[1])))
	}
}
cnvs_u=cnvs_u[-1,]
cnvs_uLocation=paste(cnvs_u$chrom, cnvs_u$start, cnvs_u$end, sep='_')
cnvs_uN=aggregate(cnvs_u, by=list(cnvs_uLocation), FUN=mean)
cnvs_uN$sec=data.frame(aggregate(cnvs_u$sec, by=list(cnvs_uLocation),function(z) paste0(sort(z), collapse=',')))[,2]
cnvs_uN$Cellular.Fraction=data.frame(aggregate(cnvs_u$Cellular.Fraction, by=list(cnvs_uLocation),function(z) paste0(sort(z), collapse=',')))[,2]
cnvs_uN$Total.Copy.Number=floor(cnvs_uN$Total.Copy.Number)
cnvs_uN$Minor.Copy.Number=round(cnvs_uN$Minor.Copy.Number)
cnvs_uN$Major.Copy.Number=cnvs_uN$Total.Copy.Number-cnvs_uN$Minor.Copy.Number
# cnvs_uN is a dataframe that lists all the unique cnvs which we believe to exist across the sections, along with the cellular fractions across those sections
cnvs_uN=cnvs_uN[order(as.numeric(cnvs_uN$chrom), cnvs_uN$start),]
cnvs_clonal=cnvs_uN[which(lapply(strsplit(cnvs_uN$Cellular.Fraction, ','),function(x) median(as.numeric(x)))<0.6),]
# cnvColDf is a dataframe that lists the CNVs across sections derived by FACETS
# }}}
# Amp-seq infered vectors of CNA
# {{{
library('stringr')
amp=read.table('~/bdata/SF10711/amp.seq/csv_files/SF10711_variant_frequencies.csv', sep=',', header=T)
amp=amp[grep('SF10711_9_', amp$sample),]
names=paste('sf', str_pad(gsub('.*_', '', amp$sample), width=3, pad=0, side='left'),sep='')
genes_to_include=c('PTPRS', 'C10orf54', 'CARD14', 'GALNT2', 'THSD7B', 'LAMA2', 'MLX', 'PLEKHA5', 'MLL4', 'KRT17', 'OR8K3', 'INA', 'RYR2', 'CDKN1A', 'ATRX', 'TP53', 'TRIB2', 'TTN', 'TMCO4', 'PPIG', 'SSH1', 'PTPRZ1', 'ZSCAN10', 'HOXD4', 'OR11L1', 'IDH1', 'LRP2')
amp=amp[,c(1, grep(paste(genes_to_include, collapse='|'), colnames(amp)))]
amp=amp[grep('SF10711_9_', amp$sample),]
colnames(amp)=gsub('C10orf54', 'VSIR', colnames(amp)) 
colnames(amp)=gsub('C14orf37', 'ARMH4', colnames(amp))
colnames(amp)=gsub('MLL4', 'KMT2B', colnames(amp))
# }}}
# CNV-kit RNA vectors of CNA
# {{{
cnvKit=fread('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/outputs/called_cnvs.csv')
cnvKit=cnvKit[,-1]
# read in imputed CNA data from all sections via cnvKit-RNA from RNA-seq
cnvKitOut=data.frame(matrix(nrow=1, ncol=ncol(cnvKit)+1)) 
closeness=0.03
j=1
for(i in seq(1,nrow(cnvKit))){
	cnvsTemp=cnvs_clonal[which(cnvs_clonal$chrom==cnvKit$chromosome[i]),]
	matchVec=((abs(cnvKit$start[i]-cnvsTemp$start)<(closeness*cnvKit$start[i])) | (abs(cnvKit$end[i]-cnvsTemp$end)<(closeness*cnvKit$end[i])))
	print(length(which(matchVec)))
	if(!(any(matchVec))){next}
	cnvKitOut[j,]=cbind(cnvsTemp$Group.1[which(matchVec)], cnvKit[i,])
	j=j+1
}
colnames(cnvKitOut)=c('Group.1', colnames(cnvKit))
cnvKitOut[,seq(5,ncol(cnvKitOut))]=t(apply(cnvKitOut[,-(1:4)], 1, function(z) smooth.spline(z, df=15)$y))
cnvKitOut=cnvKitOut[-4,]
# code below for visual inspection of whether patterns from CNV-kit match patterns from the exome data
# for(i in seq(1,nrow(cnvKitOut))){
# 	print(data.frame(cnvKitOut[i,],cnvs_clonal[which(cnvs_clonal$Group.1==cnvKitOut$Group.1[i]),]))
# }
cnvKitOut=data.frame(name=c('chr2pq_del', 'chr2q_del', 'chr11q_amp', 'chr17q_amp', 'chr18p_del'), cnvKitOut)
rownames(cnvKitOut)=cnvKitOut$name
cnvKitOut=cnvKitOut[,-seq(1,5)]
cnvKitOut=t(cnvKitOut)
cnvKitOut=cnvKitOut[which(rownames(cnvKitOut) %in% amp$sample),]
cnvKitOut=abs(cnvKitOut)
# }}}
# clonal abundance plot
# {{{
setwd('~/bdata/SF10711/clonality/figures')
library('stringr')
# first need vectors of clonal abundance from pyclone for exome and ampseq 
cluster=data.frame(fread('/mnt/bdata/patrick/SF10711/clonality/pyclone/outputs/tables/cluster.tsv'))
cluster=reshape(cluster[,c(1,2,4)], idvar='sample_id', timevar='cluster_id', direction='wide')
rownames(cluster)=cluster[,1]
cluster=cluster[,-1]
cluster[is.na(cluster)]=.99
amp=read.table('~/bdata/SF10711/amp.seq/csv_files/SF10711_variant_frequencies.csv', sep=',', header=T)
amp=amp[-(nrow(amp)),]
amp$sample=paste('sf', str_pad(gsub('.*_', '',amp$sample ), width=3, pad=0, side='left'),sep='')
cluster=cluster[order(rownames(cluster)),]
amp=amp[match(rownames(cluster), amp$sample),]
scaleFac=(amp$TP53_G338T_freq)/cluster[,3]
cluster=data.frame(apply(cluster,2, function(x) (x*scaleFac)))
clustM=apply(cluster,2,mean)
cluster=cluster[,order(clustM)]
colnames(cluster)=paste('mean', seq(4,0), sep='.')
barOut=data.frame(matrix(nrow=nrow(cluster), ncol=1))
barOut[,c(1,2,3)]=cluster[,c(1,2,3)]
barOut[,4]=cluster$mean.1-cluster$mean.2
barOut[,5]=cluster$mean.0-cluster$mean.1-cluster$mean.3-cluster$mean.4
colnames(barOut)=c('Clone 0', 'Clone 1', 'Clone 2', 'Clone 3', 'Clone 4')
barOut[barOut<0]=0
#barOut$'Clone 1'=barOut$'Clone 1'+barOut$'Clone 2'
#barOut=barOut[,-4]
#colnames(barOut)=c('Clone 2', 'Clone 4', 'Clone 5', 'Clone 3', 'Clone 1', 'Clone 0')
barOutS=apply(barOut,1, sum)
amp$TP53_G338T_freq[11]=0.73
scaleFac=(amp$TP53_G338T_freq)/barOutS
barOut=data.frame(apply(barOut,2, function(x) (x*scaleFac)))
barOut=data.frame(section=amp$sample, barOut)
barOut$section=as.factor(gsub('sf0*', '', barOut$section))
barOut=barOut[,c(1, order(colnames(barOut)[-1])+1)]
barOutP=melt(barOut)
barOutP$section<- factor(barOutP$section, levels=unique(barOutP$section))
pdf('clonal_abundance.pdf')
print(ggplot(barOutP, aes(x=section, y=value, fill=variable)) +
	geom_bar(stat="identity", width=1)+
	scale_y_continuous(expand=c(0,0))+
	labs(fill='Clone',title = "SF10711 Clonal Abundance", x = "Section", y = "Clonal frequency")+
	scale_fill_manual(values=c("#b15928", "#e31a1c", "#33a02c", "#b2df8a", "#1f78b4"), name="Clones",labels=c('Clone 0', 'Clone 1', 'Clone 2', 'Clone 3', 'Clone 4'))+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5,face = "bold", size = 20),plot.subtitle=element_text(hjust = 0.5, size=16), axis.text.x =element_text(size=6,angle = 90, vjust = 0.5, hjust=1, color='white'), axis.text.y =element_text(size=20), axis.title.x =element_text(size=20), axis.title.y =element_text(size=20),plot.margin=unit(c(.2,.5,.2,.2),"cm")))
dev.off()
# }}}
