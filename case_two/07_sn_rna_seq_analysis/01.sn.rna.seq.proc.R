# custom function to generate wide format matrix and split reads aligned to multiple sites
# {{{
library('data.table')
library('plyr')
i=1
for(ea in list.files()[grep('_R2_001.fastq.gz.extract.trim.alignAligned.sortedByCoord.out.bam.featureCounts.bam.sort.counts.tsv.gz', list.files())]){
	cat(paste('reading in ', ea, '...', sep=""))
	x=fread(ea, sep='\t', header=T)
	cat(paste('done', '\n', sep=""))
	x=dcast.data.table(x, gene~cell, value.var='count')
	x=as.data.frame(x)
	xSplit=strsplit(x$gene,',')
	xSplit1=list(lapply(xSplit, function(x) x[1]))
	xSplit1[length(xSplit1)+1]=NA
	xSplit1=unlist(xSplit1)
	xSplit2=unlist(lapply(xSplit, function(x) x[2]))
	xSplit3=unlist(lapply(xSplit, function(x) x[3]))
	genecounts=apply(x[which(unlist(lapply(xSplit, length))==1), -1], 1, sum, na.rm=T)
	names(genecounts)=x$gene[which(unlist(lapply(xSplit, length))==1)]

	geneOut=xSplit1
	threeOut=intersect(which(!(xSplit2 %in% xSplit1)), which(!(xSplit3 %in% xSplit1)))
	twoOut=setdiff(which(!(xSplit2 %in% xSplit1)), threeOut)
	threeOutalt=setdiff(which(!(xSplit3 %in% xSplit1)), union(twoOut, threeOut))
	allAlt=union(union(threeOut, twoOut), threeOutalt)
	work=which(unlist(lapply(xSplit, length))<4 & unlist(lapply(xSplit, length))>1)
	fracs=setdiff(work, allAlt)
	fracsTwo=seq(1,length(xSplit))[fracs][which(unlist(lapply(xSplit, length))[fracs]==2)]
	fracsThree=seq(1,length(xSplit))[fracs][which(unlist(lapply(xSplit, length))[fracs]==3)]

	fracsTwoAdd=data.frame(x[fracsTwo, 1], x[fracsTwo,-1]/2)
	fracsThreeAdd=data.frame(x[fracsThree, 1], x[fracsThree,-1]/3)
	fracsThreeAdd=rbind(fracsThreeAdd, data.frame(x[fracsThree, 1], x[fracsThree,-1]/3))
	
	fracsTwoAdd[,1]=xSplit2[fracsTwo]
	fracsThreeAdd[,1]=c(xSplit2[fracsThree], xSplit3[fracsThree])
	x$gene[fracsTwo]=xSplit1[fracsTwo]
	x$gene[fracsThree]=xSplit1[fracsThree]
	x[fracsTwo,-1]=x[fracsTwo,-1]/2
	x[fracsThree,-1]=x[fracsThree,-1]/3
	
	x$gene[twoOut]=xSplit2[twoOut]
	x$gene[c(threeOut, threeOutalt)]=xSplit3[c(threeOut, threeOutalt)]
	
	fracsTwoC=data.frame(matrix(nrow=length(fracsTwo), ncol=2))
	fracsTwoC[,1]=genecounts[xSplit1[fracsTwo]]
	fracsTwoC[,2]=genecounts[xSplit2[fracsTwo]]

	xSplit=strsplit(x$gene,',')
	table(unlist(lapply(xSplit, length)))
	x=x[-which(unlist(lapply(xSplit, length))>1),]
	colnames(x)=paste(colnames(x), gsub("_.+", "", ea), sep=".")
	x[,1]=sapply(x[,1], function(x) gsub("\\..+", "", x))
	colnames(x)[1]="Gene"
	x=data.table(x)
	x=x[,lapply(.SD, sum, na.rm=T), by=list(Gene)]
	x=data.frame(x)
	x[,1]=as.factor(x[,1])
	if(i==1){
		expr.raw=x
		genes=expr.raw[,1]
	} else {
		genes=union(expr.raw[,1], x[,1])
		expr.raw=expr.raw[match(genes, expr.raw[,1]),]
		expr.raw[,1]=genes
		expr.raw=cbind(expr.raw, x[match(genes, x[,1]),-1])
	}
	i=i+1
}
expr.raw[is.na(expr.raw)]=0
read.per.cell=apply(expr.raw[,-1], 2, sum)
for (ea in seq(1,12)){
	print(ea)
	print(sum(read.per.cell[grep(paste("N7.*", ea, sep=""), names(read.per.cell))])/1E5)
}
read.per.cell=apply(expr.raw[grep('ERCC', expr.raw$Gene),-1], 2, sum)
ercc_check=data.frame(matrix(nrow=1, ncol=13))
i=1
for(i_row in grep('ERCC', expr.raw$Gene)){
	for (ea in seq(1,12)){
		temp=expr.raw[i_row,]
		ercc_check[i, ea]=sum(temp[grep(paste("N7.*", ea, sep=""), names(temp))])
	}
	print(i_row)
	i=i+1
}
gene_check=data.frame(matrix(nrow=1, ncol=13))
i=1
for(i_row in seq(1,1000)){
	for (ea in seq(1,12)){
		temp=expr.raw[i_row,]
		gene_check[i, ea]=sum(temp[grep(paste("N7.*", ea, sep=""), names(temp))])
	}
	print(i_row)
	i=i+1
}
write.csv(expr.raw, file="expr.ensg.counts.csv", row.names=F)
# }}}
# tidying expression matrix
# {{{
library('data.table')
gtf=fread('/home/shared/hg_align_db/./GRCh38.p13.2019-02-28.gencode/gencode.v33.chr_patch_hapl_scaff.annotation.ERCC.gtf', skip=5)
colnames(gtf)=c("chr", "db", "type", "start", "end", "score", "strand", "frame", "extra")
gtf=gtf[which(gtf$type=="gene"),]
gtf.ENSG=gsub(".*gene_id\\s\"\\\"*|\\..*", "", gtf$extra)
gtf.type=gsub(".*gene_type\\s\"\\\"*|\".*", "", gtf$extra)
gtf.name=gsub(".*gene_name\\s\"\\\"*|\".*", "", gtf$extra)
gtf.out=data.frame(gtf[,1], gtf.ENSG, gtf.type, gtf.name)
expr=fread("expr.ensg.counts.csv")
expr=data.table(gtf.out[match(expr$Gene, gtf.out$gtf.ENSG),c(3,4)], expr)
expr$gtf.type=as.character(expr$gtf.type)
expr$gtf.name=as.character(expr$gtf.name)
expr$gtf.type[which(is.na(expr$gtf.type))]="NA"
expr$gtf.name[which(is.na(expr$gtf.name))]="NA"
expr$gtf.type[grep('^ERCC-', expr$Gene)]="ERCC"
expr$gtf.name[grep('^ERCC-', expr$Gene)]=expr$Gene[grep('^ERCC-', expr$Gene)]
expr=expr[, lapply(.SD, sum), by = gtf.name, .SDcols = !c("Gene", "gtf.type")]
expr=data.table(gtf.out$gtf.type[match(expr$gtf.name, gtf.out$gtf.name)], expr)
detec.p.gene=apply(expr[,-c(1,2)], 1, function(x) length(which(x>0)))
# expr=expr[which(detec.p.gene>80),]
sif=data.frame(colnames(expr)[-c(1,2)], substr(colnames(expr)[-c(1,2)],1,12), substr(colnames(expr)[-c(1,2)],14,20))
colnames(sif)=c('name', 'cell.barc', 'plate.barc')
sif=data.frame(sif, all_Patients=rep("all_cells", nrow(sif)))
# N701,N702,etc.
indeces=c("TAAGGCGA", "CGTACTAG", "AGGCAGAA", "TCCTGAGC", "GGACTCCT", "TAGGCATG", "CTCTCTAC", "CAGAGAGG", "GCTACGCT", "CGAGGCTG", "AAGAGGCA", "GTAGAGGA")
index.order=c(1, 2, 9, 12, 3,4, 5,6, 10, 7,8,11)
indeces=indeces[index.order]
spike.ind=expr[grep("ERCC", expr$gtf.name),1]
sif=sif[order(match(sif$plate.barc, indeces)),]
sif[,5]=rep("NA", nrow(sif))
sif[grep("N701|N702|N709", sif$plate.barc),5]=17
sif[grep("N712|N703|N704", sif$plate.barc), 5]=53
sif[grep("N705|N706|N710", sif$plate.barc), 5]=93
sif[grep("N707|N708|N711", sif$plate.barc), 5]=117
colnames(sif)[5]="slice"
sif=sif[order(as.numeric(as.character(sif$slice)), as.numeric(sif$plate.barc)),]
sif[,5]=as.factor(sif[,5])
expr=data.frame(expr)
expr=expr[,c(1,2,match(sif$name, colnames(expr)))]
sif$name=paste(sif$plate.barc, sif$cell.barc, sep=".")
write.csv(sif, file='sif.csv', row.names=F)
colnames(expr)[-c(1,2)]=sif$name
colnames(expr)[1:2]=c("class", "Gene")
write.csv(expr, file="expr.ensg.name.counts.csv", row.names=F)
# }}}
# quality control plots
# {{{
library('RUVSeq')
library('RColorBrewer')
library('future')
library('WGCNA')
cl=makeCluster(detectCores()-5)
plotKadj=function(dfin, k, randCells){
	dfinBak=dfin
	genes.sum=parApply(cl=cl, dfin, 1, sum)
	dfin=dfin[order(genes.sum, decreasing=T)[1:10000],]
	addition=min(parApply(cl=cl, dfin,2,min))
	if(addition<0){addition=abs(addition)}else{addition=0}
	exprPR=prcomp(t(log2(dfin+1+addition)))
	exprPRplot=data.frame(exprPR$x, color=gsub('\\..*', '', rownames(exprPR$x)))
	distO=cor(t(dfin))
	samps=gsub('\\..*', '', colnames(dfin))
	aovPSamps=summary(aov(exprPR$x[,1] ~ as.factor(samps)))
	aovPPval=data.frame(aovPSamps[[1]])[1,5]
	aovPMalig=summary(aov(exprPR$x[,1] ~ as.factor(amps)))
	aovMPval=data.frame(aovPMalig[[1]])[1,5]
	exprRLE=parApply(cl=cl, dfin, 2, function(x) log10(x/(median(x,na.rm=T)+1)+1))
	exprRLEplot=melt(data.frame(exprRLE[,randCells]))
	exprRLEplot[,2]=as.numeric(exprRLEplot[,2])
	distO=distO[upper.tri(distO)]
	readsSum=parApply(cl=cl, dfinBak[,randCells], 2, sum)
	readsSum=readsSum/max(readsSum)*100
	pa=ggplot(exprRLEplot, aes(x=variable, y=value))+geom_boxplot(outlier.shape=NA, notch=F)+labs(title=paste('RUVg RLE plot, k=', k), subtitle=paste('Median:', signif(median(unlist(dfin)), 4), ', Median variance: ', signif(median(apply(dfin, 2, var)),4)), x='Nuclei', y='RLE (log count/median)')+theme_classic()+theme(legend.position="none", axis.text.x=element_blank())
	pb=ggplot(data.frame(value=distO), aes(x=value))+geom_histogram(color='black', fill='white',bins=30)+labs(y='Counts', x='Pearson correlation', title='Distribution of Pearson correlations', subtitle=paste('ANOVA Barcode Pval:' , signif(aovPPval,3), '\nANOVA Malignancy Pval:', signif(aovMPval,3))) + theme_classic()+ theme(legend.position="none" )
	pc=ggplot(data.frame(barcodes=names(readsSum), reads=readsSum), aes(x=barcodes, y=reads, fill=barcodes))+geom_bar(stat='identity')+labs(x='Barcodes', y='Relative percentage read depth for sample')+theme_classic() + theme(legend.position="none", axis.text.x=element_blank())
	pd=ggplot(exprPRplot, aes(x=PC1, y=PC2, color=color))+geom_point(show.legend=F)+labs(title=paste('PCA plot, PC1 v PC2, k=', k))+theme_classic()
	pe=ggplot(exprPRplot, aes(x=PC1, y=PC3, color=color))+geom_point(show.legend=F)+labs(title=paste('PCA plot, PC1 v PC3, k=', k))+theme_classic()
	pf=ggplot(exprPRplot, aes(x=PC2, y=PC3, color=color))+geom_point(show.legend=F)+labs(title=paste('PCA plot, PC2 v PC3, k=', k))+theme_classic()
	multiplot(pa,pc,pb,pd,pe,pf, cols=2)
}
# }}}
# perform RUVg adjustment based on ERCC controls
# {{{
expr.n=data.frame(fread("expr.ensg.name.counts.csv"))
expr.n.b=expr.n[,c(1,2)]
rownames(expr.n)=expr.n$Gene
spikes=expr.n$Gene[grep("^ERCC-", expr.n$Gene)]
randCells=sample(seq(1,ncol(expr.n)-2), 70)
pdf('samples_adjustment_dashboards_alt.pdf')
plotKadj(dfin=expr.n[,-c(1,2)], k='raw', randCells=randCells)
for(kn in c(seq(0,10,2),seq(15, 25, 5), seq(30, 100,10))){
	if(kn==0){
		 plotKadj(expr.n[,-c(1,2)], kn, randCells=randCells)
	}
	else{
		exprRun=RUVg(as.matrix(round(expr.n[,-c(1,2)])), spikes, k=kn, round=F, center=TRUE, isLog=F)
		plotKadj(exprRun$normalizedCount, kn, randCells=randCells)
	}
	print(kn)
}
expr.n=RUVg(as.matrix(round(expr.n[,-c(1,2)])), spikes, k=10, round=F, center=TRUE, isLog=F)$normalizedCounts
expr.n=data.frame(expr.n.b, expr.n)
expr.n=expr.n[-grep('^ERCC-', expr.n$Gene),]
expr.n.b=expr.n[,c(1,2)]
expr.n=expr.n[,-c(1,2)]
addition=min(parApply(cl=cl, expr.n,2,min))
if(addition<0){addition=abs(addition)}else{addition=0}
expr.n=expr.n+addition
write.csv(data.frame(expr.n.b, expr.n), file="expr.ensg.name.ruvg.k10.unrounded.csv", row.names=F)
reads.sum=parApply(cl=cl, expr.n, 2, sum)
for(ea in seq(1,ncol(expr.n))){
	expr.n[,ea]=expr.n[,ea]/(reads.sum[ea]/1E6)
}
plotKadj(dfin=expr.n[,], k='cpm', randCells=randCells)
expr.n=log2(expr.n+1)
expr.n=data.frame(expr.n.b, expr.n)
dev.off()
write.csv(expr.n, file='expr.ensg.name.ruvg.k10.unrounded.cpm.log2.csv', row.names=F)
################################################################################
