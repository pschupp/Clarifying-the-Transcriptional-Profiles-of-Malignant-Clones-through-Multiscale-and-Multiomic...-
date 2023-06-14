facets_run=function(WD, cval1, sub=NA){
	library('facets')
	library('future')
	plan(multicore)
	set.seed(666)
	setwd(WD)
	out=list()
	k=1
	if(is.na(sub))	{itFiles=seq(1, length(list.files('pileup')))} else
					{itFiles=grep(sub, list.files('pileup'))}
	for(i in itFiles){
		out[[k]] <- future({
			out_frame=data.frame(matrix(nrow=1, ncol=5))
			j=1
			ea=list.files('pileup')[i]
			snp <- readSnpMatrix(paste0('pileup/',ea))
			xx <- preProcSample(snp, ndepth=10, het.thresh=0.1, ndepthmax=10000, cval=cval1)
			pdf(paste0('facets/',gsub('\\..*', '', ea),'.pdf'))
			for(cval in seq(floor((cval1/100)+1)*100, 500, 50)){
				name=paste(gsub('\\..*', '', ea), 'cval', cval, sep='.')
				oo <- procSample(xx, cval = cval, dipLogR=0)
				fit <- emcncf(oo)
				out_frame[j,]=data.frame(cval=cval, loglike=fit$loglik, purity=fit$purity, ploidy=fit$ploidy, dipLogR=fit$dipLogR)
				j=j+1
				plotSample(oo, fit, sname = name, plot.type='em')
				fit$cncf$`Major Copy Number` <- fit$cncf$tcn.em - fit$cncf$lcn.em
				colnames(fit$cncf)[is.element(colnames(fit$cncf), c("cf.em", "tcn.em", "lcn.em"))] <- c("Cellular Fraction", "Total Copy Number", "Minor Copy Number")
				fwrite(fit$cncf, file = paste0('facets/',name,'.csv'))
			}
			dev.off()
			colnames(out_frame)=c('cval', 'logLike', 'purity', 'ploidy', 'dipLogR')
			fwrite(out_frame, file=paste0('facets/',gsub('\\..*', '', ea) ,'.stats.csv'))
		})
		k=k+1
	}
}
cnAssign=function(cnv,name){
	out=data.frame(matrix(nrow=1,ncol=3))
	colnames(out)=c('Cellular.Fraction', 'Major.Copy.Number', 'Minor.Copy.Number')
	for(i in seq(1, nrow(gtf))){
		temp=(cnv[intersect(intersect(which(cnv$chrom==gtf$chr[i]),which(cnv$start<gtf$mean[i])),which(cnv$end>gtf$mean[i])),])
		if(length(temp$Cellular.Fraction)>0){
			out=rbind(out,temp[,c(12,15,14)])
		} else{
			out=rbind(out, rep(NA,3))
		}
	}
	out=out[-1,]
	colnames(out)=paste(name, colnames(out), sep='.')
	gtf=cbind(gtf,out)
	gtf
}
citupToTree=function(citDir, snpDir,outDir){
	library('treemap')
	library('data.tree')
	library('ape')
	library('ggtree')
	library('DiagrammeRsvg')
	library('DiagrammeR')
	setwd(citDir)
	objValue=c()
	for(ea in list.files('tmp/tmp/tree/')){
		setwd(citDir)
		setwd(paste0('tmp/tmp/tree/', ea))
		if(length(grep('results', list.files()))==0){next}
		#print(ea)
		objValue=c(objValue, as.numeric(readLines('results', n=2)[2]))
	}
	setwd(snpDir)
	if(length(grep('cluster', list.files()))==1){
		loci=fread('loci.tsv', sep='\t', header=T)
		cluster=fread('cluster.tsv', sep='\t', header=T)
	} else if(length(grep('.*output.tsv', list.files()))==1){
		output=fread(list.files()[grep('.*output.tsv', list.files())])
		output$size=table(output$cluster_id)[match(output$cluster_id, data.frame(table(output$cluster_id))$Var1)]/length(unique(output$sample_id))
		cluster=output[,c(2, 3, 7, 4,5)]
		colnames(cluster)[4:5]=c('mean','std')
		cluster=aggregate(cluster, by=list(paste0(cluster$sample_id, cluster$cluster_id)), mean)
		cluster$sample_id=cluster$Group.1
		cluster=cluster[,-1]
		loci=output[,c(seq(1,5))]
		lociNew=data.frame(matrix(ncol=ncol(loci), nrow=1))
		colnames(lociNew)=colnames(loci)
		for(clo in unique(loci$cluster_id)){
			temp=loci[which(loci$cluster_id==clo),]
			temp=temp[order(temp$cellular_prevalence, decreasing=T)[sample(seq(1,nrow(temp)),10)],]
			lociNew=rbind(lociNew,temp)
		}
		lociNew=lociNew[-1,]
		loci=lociNew
	}
	x=aggregate(loci$cellular_prevalence, by=list(loci$cluster_id), FUN=mean)
	x$ord=rank(x$x)
	x$ord=rank(1-x$x)
	loci$cluster_id=x$ord[match(loci$cluster_id, x$Group.1)]
	loci=loci[order(loci$cluster_id,as.numeric(gsub(':.*','', loci$mutation_id)), as.numeric(gsub('.*:','', loci$mutation_id))),]
	setwd(citDir)
	l=1
	library('future')
	plan(multicore, workers=10)
	p=1
	out=list()
	loci$cluster_id=loci$cluster_id-2
	for(ea in list.files('tmp/tmp/tree/')){
		setwd(citDir)
		eaOg=ea
		setwd(paste0('tmp/tmp/tree/', ea))
		if(length(grep('results', list.files()))==0){next}
		print(ea)
		x=readLines('results', n=15)[-seq(1,3)]
		score=round(as.numeric(readLines('results', n=2)[2]),3)
		x=x[-seq(grep('gamma_matrix', x), length(x))]
		x=as.data.frame(t(sapply(strsplit(x, ' '), unlist)))
		pcTable=data.frame(matrix(nrow=length(unique(unlist(x))), ncol=3))
		i=1
		for(eaX in unique(unlist(x))){
			parent=paste(x[which(x[,2]==eaX),1], collapse=',')
			child=paste(x[which(x[,1]==eaX),2], collapse=',')
			if(length(parent)==0){parent=NA}
			if(length(child)==0){child=NA}
			pcTable[i,]=c(eaX,parent,child)
			i=i+1
		}
		colnames(pcTable)=c('clone', 'parent', 'child')
		parent=pcTable$clone[which(nchar(pcTable$parent)==0)]
		pPrime=99
		outF=c()
		for(tail in pcTable$clone[which(nchar(pcTable$child)==0)]){
			pPrime=i=99
			out=c()
			while(pPrime!=parent){
				pPrime=pcTable$parent[which(pcTable$clone==tail)]
				if(i==99){out=c(pPrime, tail, out)}else {out=c(pPrime, out)}
				tail=i=pPrime
			}
			outF=c(outF, paste(out, collapse='/'))
			#print(outF)
		}
		outF=data.frame(pathString=outF)
		outF=apply(outF,1,function(x) paste0('-1/', x))
		outF=data.frame(pathString=outF)
#		if(nrow(outF)==1){next}
		outNode=as.Node(outF)
		labelVec=data.frame(matrix(nrow=7, ncol=1))	
		bigIter=as.numeric(unique(unlist(strsplit(as.character(paste(unlist(outF), collapse='/')), '/'))))
		for(ea in bigIter[-length(bigIter)]){
			temp=loci[which(loci$cluster_id==ea),]
			temp=temp[which(temp$cellular_prevalence>0),]
			temp2=aggregate(temp[,c('cluster_id', 'cellular_prevalence')], by=list(temp$mutation_id), FUN=mean)
			temp2=temp2[order(temp2$cellular_prevalence, decreasing=T),]
			labelVec[[(ea+2)]]=temp2$Group.1[1:7]
		}
		labelVecL=lapply(labelVec, paste, collapse='\n')
		labelVecL=lapply(labelVecL, function(x) gsub('\nNA', '', x))
		outNode$Do( function(node) SetEdgeStyle(node, inherit=FALSE, arrowhead = "vee", color = "grey35", penwidth = 2, label=labelVecL[(as.numeric(node$name)+1)]))
		setwd(outDir)
		treeAsSVG <- export_svg(render_graph(ToDiagrammeRGraph(outNode),title=paste0('Score:', score, ',', eaOg)))
		writeLines(treeAsSVG, paste0('file', eaOg, '.svg'))
		l=l+1
#	})
	l=l+1
#	p=p+1
	}
	setwd(outDir)
	system('for ea in *svg; do sem -j 15 inkscape $ea --export-pdf=$ea.pdf; done')
	system('qpdf --empty --pages *.pdf -- out.pdf')
}
cnAssignNumr=function(mut, cnv){
	major_cn=minor_cn=rep(NA, nrow(mut))
	for(i in seq(1,nrow(cnv))){
		ind=which(mut$Chromosome==cnv$chrom[i] & mut$Position>=cnv$start[i] & mut$Position<=cnv$end[i])
		if(length(ind)==0){next}
		major_cn[ind]=rep(cnv$'Major Copy Number'[i], length(ind))
		minor_cn[ind]=rep(cnv$'Minor Copy Number'[i], length(ind))
	}
	major_cn[is.na(major_cn)]=1
	minor_cn[is.na(minor_cn)]=1
	mut=data.table(mut, major_cn=major_cn, minor_cn=minor_cn)
	return(mut)
}
