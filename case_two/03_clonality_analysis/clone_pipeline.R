# Goal is to create a concencus tree of mutations with cluster prevalance tracked through the ~85 sections, also compare to primary tumor. For primary tumor goal is to show that the derived clones are closer to trunk of tree
#{{{
# GENERATING MPILEUP FILES
cd /mnt/bdata/patrick/SF10711/exome.seq/
nice -n -18 samtools mpileup -B -q 1 -f /home/shared/hg_align_db/GRCh37.38/GRCh37.primary_assembly.genome.fa bamfiles/Z00346.bwa.realigned.rmDups.recal.bam | nice -n -18 varscan mpileup2snp --min-freq-for-hom 0.95 --min-var-freq 0.05 --p-value 0.1 --output-vcf 1  > vcf_file/Z00346.bwa.realigned.rmDups.recal.bam.vcf
# Only SNPs will be reported
# Min coverage:   8
# Min reads2:     2
# Min var freq:   0.05
# Min avg qual:   15
# P-value thresh: 0.1
# Reading input from Z00346.bwa.realigned.rmDups.recal.bam.pileup
# 166556 Z00346.bwa.realigned.rmDups.recal.bam.vcf
# 65876 Blood_kapa-NGv3-RPE100-NGv3-gatk-haplotype.vcf
cd /mnt/bdata/patrick/SF10711/exome_seq/
for ea in bamfiles/SF10711_9*bam; do
	echo $ea
	name=$(sed 's/.*\///g' <<< $ea)
	sem -j 4 nice -n -18 snp-pileup --count-orphans --max-depth=10000 --gzip --progress --min-map-quality=1  --min-base-quality=20 --min-read-counts 8,0 vcf_file/Z00346.bwa.realigned.rmDups.recal.bam.vcf pileup/$name.pileup bamfiles/Z00346.bwa.realigned.rmDups.recal.bam $ea
done

cd /mnt/bdata/patrick/SF10711/exome_seq/
for ea in bamfiles/Z00363*bam; do
	echo $ea
	name=$(sed 's/.*\///g' <<< $ea)
	sem -j 4 nice -n -18 snp-pileup --count-orphans --max-depth=10000 --gzip --progress --min-map-quality=1  --min-base-quality=20 --min-read-counts 8,0 vcf_file/Z00346.bwa.realigned.rmDups.recal.bam.vcf pileup/$name.pileup bamfiles/Z00346.bwa.realigned.rmDups.recal.bam $ea
done
#  Calculating SNP count...done.
# Finished in 598.168231 seconds.========================================] 100 %
# double free or corruption (out)
# Calculating SNP count...done.
# Finished in 621.123145 seconds.========================================] 100 %
# double free or corruption (out)
# Calculating SNP count...done.
# Finished in 649.705248 seconds.========================================] 100 %
# double free or corruption (out)
# Calculating SNP count...done.
# Finished in 652.782948 seconds.========================================] 100 %
# double free or corruption (out)
#}}}
# running FACETS
#{{{
source('/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R')
#facets_run("/mnt/bdata/patrick/SF10711/exome.seq/", 50)
facets_run("/mnt/bdata/patrick/SF10711/exome_seq/", 150, sub='SF10711')
#}}}
# convert pilup locations to genes and protein changes
#{{{
cd /mnt/bdata/patrick/SF10711/exome_seq/pileup/
for(ea in list.files()[grep('Z00363.*gz', list.files())]){
	pileup=fread(ea)
	pileup=pileup[-which(pileup$Chromosome %in% c('chrM', 'chrY','chrX')),]
	pileupDepth=apply(pileup[,c('File2R','File2A')], 1, sum)
	pileup=pileup[(pileupDepth>200),]
	pileOut=data.frame(chr=gsub('chr', '', pileup$Chromosome), start=pileup$Position, stop=pileup$Position, allele=paste(pileup$Ref, pileup$Alt, sep='/'))
	fwrite(pileOut, row.names=F, col.names=F, sep='\t', quote=F, file=paste0(gsub('gz', '', ea), 'vep.gz'), compress='gzip')
}
for ea in Z00363*vep.gz; do 
	eaN=${ea/.*/}
	/opt/vep_ensembl/ensembl-vep/vep	-i $ea -o $eaN.vep.txt --cache --assembly GRCh37 --offline --dir_cache /home/shared/vep_cache/ \
										--fasta /home/shared/vep_cache/fasta/GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa --hgvs --hgvsg --shift_hgvs 1 \
										--symbol --refseq --fork 1 --pick --tab \
										--fields  "Uploaded_variation,SYMBOL,STRAND,Feature,HGVSc,HGVSp,HGVSg" 
done
#}}}
# CNV-kit RNA workflow
# FOR PYCLONE, USING AMPSEQ SNPS + CNVkit RNA-seq CNVS CALIBRATED BY EXOME CNVS
# MAKE DATAFRAME OF CALLED CNVS FROM EXOME-SEQ+FACETS
#{{{
setwd('/mnt/bdata/patrick/SF10711/cnv.analysis/rna.seq/cnvkit-rna')
# read in exome cnv data from FACETS
cnv=list()
cnv[[1]]=fread('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-1-22.cval.450.csv')
cnv[[2]]=fread('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-1-46.cval.450.csv')
cnv[[3]]=fread('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-2-85.cval.500.csv')
cnv[[4]]=fread('/mnt/bdata/patrick/SF10711/exome_seq/facets/SF10711_9-2-123.cval.450.csv')
names(cnv)=c('s22', 's46', 's85', 's123')
# get concensus cnvs
cnv=lapply(cnv, function(x) x=x[which(x$chrom %in% seq(1,22)),])
cnv=lapply(cnv, function(x) x=x[which(x$'Minor Copy Number'!=1|x$'Major Copy Number'!=1),])
cnv=lapply(cnv, function(x) x=x[which((x$end-x$start)>1E7),])
cnv=lapply(cnv, function(x) x=x[,c(1,6,10,11,12,13,14,15)])
for(i in seq(1, length(cnv))){
	cnv[[i]]$section=rep(names(cnv)[i], nrow(cnv[[i]]))
}
cnvs=as.data.frame(data.table::rbindlist(cnv))
cnvs=cnvs[order(as.numeric(cnvs$chrom), cnvs$start, cnvs$section),]
cnvs=rbind(cnvs, rep(0, ncol(cnvs)))
cnvs$clone=rep(NA, nrow(cnvs))
j=1
for(i in seq(1,(nrow(cnvs)-1))){
	x=1
	if(cnvs$chrom[i]==cnvs$chrom[i+1] & cnvs$end[i]>=(cnvs$end[i+1]*.9) & cnvs$end[i]<=(cnvs$end[i+1]*1.1)){
		cnvs$clone[i]=j
		x=0
		next
	}
	cnvs$clone[i]=j
	j=j+1
}
cnvs=cnvs[-which(as.numeric(cnvs$clone) %in% which(table(cnvs$clone)==1)),]
cnvs=cnvs[-which(is.na(as.character(cnvs$clone))),]
write.csv(cnvs, file='/mnt/bdata/patrick/SF10711/exome_seq/facets_called_cnvs.csv')
# manual compilation of CNVs
#}}}
# AGGREGATE AND GIVE NAMES TO CNVS
#{{{
cnvs=fread('/mnt/bdata/patrick/SF10711/exome_seq/facets_called_cnvs.csv')
for(cNu in cnvs$clone){
	cInd=which(cnvs$clone==cNu)
	cnvs$start[cInd]=mean(cnvs$start[cInd])
	cnvs$end[cInd]=mean(cnvs$end[cInd])
}
cnvsAg=cbind(aggregate(cnvs[,c(1,3,4,6,7,8)], by=list(cnvs$clone), function(x) round(mean(x)))[,-c(1,2)], reshape(cnvs[,c(1,5,9,10)], timevar='section', idvar=c('clone','chrom'), direction='wide'), aggregate(cnvs[,5], by=list(cnvs$clone), FUN=mean))
cnvsAg=cnvsAg[,    c(6,       1,     2,      7,       13,      3,     5,        4,      10,   9,      11,    8)]
colnames(cnvsAg)=c('chr', 'start', 'end', 'clone', 'avg_cf', 'TCN', 'MajCN', 'MinCN', 's22', 's46', 's85', 's123')
cnvsAg[,grep('^s[0-9]', colnames(cnvsAg))][is.na(cnvsAg[,grep('^s[0-9]', colnames(cnvsAg))])]=0.0
cnvsAg$event=rep(NA, nrow(cnvsAg))
cnvsAg$event[which(cnvsAg$TCN<2)]=paste0(abs(2-cnvsAg$TCN[which(cnvsAg$TCN<2)]), ' x del')
cnvsAg$event[which(cnvsAg$TCN>2)]=paste0(abs(2-cnvsAg$TCN[which(cnvsAg$TCN>2)]), ' x amp')
cnvsAg$event[which(cnvsAg$TCN==2)]='CN LOH'
#}}}
# read in gtf file which will allow for conversion between gene name and ensembl gene id
#{{{
gtf=fread('/home/shared/hg_align_db/GRCh37.38/gencode.v38lift37.annotation.gene_only.gtf')
gtf=gtf[,c(1,4,5,7,9)]
gtf=gtf[-grep('chrM|chrX|chrY', unlist(gtf[,1])),]
gtf$ensembl=gsub('\\..*', '', unlist(gsub('gene_id\ "', '', unlist(gtf[,5]))))
gtf$hugo=gsub('".*', '', unlist(gsub('.*gene_name\ "', '', unlist(gtf[,5]))))
gtf$type=gsub('".*', '', unlist(gsub('.*gene_type\ "', '', unlist(gtf[,5]))))
gtf=gtf[,-5]
colnames(gtf)=c('chr', 'start', 'end', 'strand', 'ensembl', 'gene', 'type')
gtf$chr=gsub('chr', '', gtf$chr)
gtf=as.data.frame(gtf)
gtf[,c(1,2,3)]=apply(gtf[,c(1,2,3)], 2, as.numeric)
# READ IN RNA SEQ DATA
rna=fread('/mnt/bdata/patrick/SF10711/rna.seq/raw.counts/RNAseq_non-normalized-hg19-ERCC92_trimmed_reads_stranded_8lanes_featurecounts_Q1_s2.csv')
rna=rna[,-6]
rna$Chr=gsub(';.*', '', rna$Chr)
rna$Chr=gsub('chr', '', rna$Chr)
rna$Start=gsub(';.*', '', rna$Start)
rna$End=gsub(';.*', '', rna$End)
rna$Strand=gsub(';.*', '', rna$Strand)
rna=rna[which(rna$Chr %in% seq(1,22)),]
colnames(rna)[1:5]=c('gene', 'chr', 'start', 'end', 'strand')
rna=data.frame(rna)
rna=rna[-which(apply(rna[,grep('SF10711', colnames(rna))], 1, var)==0),]
rna[,c(2,3,4,seq(6, ncol(rna)))]=apply(rna[,c(2,3,4,seq(6, ncol(rna)))], 2, as.numeric)
rna=rna[which(apply(rna[,-seq(1,5)], 1, function(x) length(which(x>10)))>94),]
# GENE COUNTS - THE GENE ENSEMBL IDS AND PER-GENE READ COUNTS CAN BE READ FROM A SIMPLE 2-COLUMN, TAB-DELIMITED FILE. 
rnaOut=data.frame(matrix(nrow=1, ncol=ncol(rna)))
# NEED TO CREATE CNV-EXPRESSION CORRELATION COEFFICIENTS – CALCULATED PER GENE VIA THE SCRIPT CNV_EXPRESSION_CORRELATE.PY INCLUDED WITH CNVKIT. 
outCor=data.frame(matrix(nrow=1, ncol=6))
colnames(outCor)=c('gene', 'chr','start', 'end', 'strand', 'Pearson.Corr')
colnames(rnaOut)=colnames(rna)
pdf('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/plots/dist_cor_genes_in_cnvs.pdf')
for(i in seq(1,nrow(cnvsAg))){
	chr=cnvsAg$chr[i]	
	start=cnvsAg$start[i]	
	end=cnvsAg$end[i]	
	rnaT=rna[which(rna$chr==chr),]
	rnaT=rnaT[intersect(which(rnaT$start>start), which(rnaT$end<end)),]
	if(nrow(rnaT)==0){
		print('skip')
		next
	}
	if(round(cnvsAg$avg_cf[i], 2)>0.6){next}
	if(nrow(rnaT)<80){next}
	distP=data.frame(vals=as.numeric(cor(t(rnaT[, grep('_22|_46|_85|_123', colnames(rnaT))]), t(cnvsAg[i, grep('s22|s46|s85|s123', colnames(cnvsAg))]))))
	print(ggplot(distP, aes(x=vals)) + geom_histogram() + theme_bw() +xlim(-1,1) + labs(title=paste0('chr', cnvsAg$chr[i],':', cnvsAg$start[i], '-', cnvsAg$end[i],', length (Mb):', (cnvsAg$end[i]-cnvsAg$start[i])/1E6), subtitle=paste0(cnvsAg$event[i], ', avg CF:', round(cnvsAg$avg_cf[i], 2), ', genes:', nrow(rnaT)), x='Pearson Correlation'))
	outCorT=data.frame(rnaT[,c(1,2,3,4,5)], Pearson.Corr=as.numeric(cor(t(rnaT[, grep('_22|_46|_85|_123', colnames(rnaT))]),t(cnvsAg[i, grep('s22|s46|s85|s123', colnames(cnvsAg))]))))
	outCor=rbind(outCor,outCorT)
	rnaOut=rbind(rnaOut, rnaT)
}
dev.off()
rnaOut=rnaOut[-1,]
outCor=outCor[-1,]
rnaOut=rna
rnaOut[,-c(1,5)]=apply(rnaOut[,-c(1,5)], 2, as.numeric)
#}}}
#  GENE INFO – A TABLE OF PER-TRANSCRIPT AND PER-GENE METADATA EXPORTED FROM ENSEMBL BIOMART. A SNAPSHOT OF THIS FILE FOR THE HUMAN TRANSCRIPTOME IS BUNDLED WITH CNVKIT UNDER DATA/ENSEMBL-GENE-INFO.HG38.TSV. [CAN PROBABLY USE THIS ONE]
#{{{
library('mygene')
info=fread('/opt/cnvkit/data/ensembl-gene-info.hg38.tsv')
conv=gtf[match(rnaOut$gene, gtf$gene),]
conv=conv[-which(is.na(conv$gene)),]
rnaOut$ENSG=conv$ensembl[match(rnaOut$gene, conv$gene)]
xli=rnaOut$gene[which(is.na(rnaOut$ENSG))]
out <- data.frame(queryMany(xli, scopes="symbol", fields="ensembl.gene", species="human", returnall=F, return.as='DataFrame'))
out$ensembl=as.character(lapply(lapply(out$ensembl, unlist), function(x) x[[1]][1]))
out=out[-intersect(which(duplicated(out$query) | duplicated(out$query, fromLast=T)), grep('^[0-9]', out$X_id)),]
out=out[match(rnaOut$gene, out$query),]
rnaOut$ENSG[which(is.na(rnaOut$ENSG))]=out$ensembl[which(is.na(rnaOut$ENSG))]
rnaOut$ENSG[which(rnaOut$ENSG=='NULL')]=NA
convFile=fread('/home/rebecca/omicon/mapping_tables/org.Hs_SYMBOL_ALIAS.csv')
rnaOut$gene[which(rnaOut$gene %in% convFile$ALIAS)]=convFile$SYMBOL[match(rnaOut$gene[which(rnaOut$gene %in% convFile$ALIAS)], convFile$ALIAS)]
xli=rnaOut$gene[which(is.na(rnaOut$ENSG))]
out <- data.frame(queryMany(xli, scopes="symbol", fields="ensembl.gene", species="human", returnall=F, return.as='DataFrame'))
out$ensembl=as.character(lapply(lapply(out$ensembl, unlist), function(x) x[[1]][1]))
out=out[-intersect(which(duplicated(out$query) | duplicated(out$query, fromLast=T)), grep('^[0-9]', out$X_id)),]
out=out[match(rnaOut$gene, out$query),]
rnaOut$ENSG[which(is.na(rnaOut$ENSG))]=out$ensembl[which(is.na(rnaOut$ENSG))]
rnaOut=rnaOut[,c(seq(1,5), 102, seq(6,101))]
rnaOut=rnaOut[!duplicated(rnaOut$ENSG),]
#}}}
# NOW READ IN NORMAL BRAIN FRONTAL LOBE DATA FROIM ALLEN BRAIN ATALS
#{{{
library('mygene')
library('pheatmap')
allen1=data.frame(fread('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/normal/allen_brain_H0351.2001/RNAseqCounts.csv'))
allen1Sif=data.frame(fread('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/normal/allen_brain_H0351.2001/SampleAnnot.csv'))
allen1Ind=intersect(which(allen1Sif$main_structure=="FL"), which(allen1Sif$replicate_sample=='No'))
allen1=allen1[,c(1,(allen1Ind+1))]
colnames(allen1)=c('gene', paste(gsub('_.*', '', allen1Sif$RNAseq_sample_name[allen1Ind]), allen1Sif$sub_structure[allen1Ind],sep='_'))
allen2=data.frame(fread('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/normal/allen_brain_H0351.2002/RNAseqCounts.csv'))
allen2Sif=data.frame(fread('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/normal/allen_brain_H0351.2002/SampleAnnot.csv'))
allen2Ind=intersect(which(allen2Sif$main_structure=="FL"), which(allen2Sif$replicate_sample=='No'))
allen2=allen2[,c(1,(allen2Ind+1))]
colnames(allen2)=c('gene', paste(gsub('_.*', '', allen2Sif$RNAseq_sample_name[allen2Ind]), allen2Sif$sub_structure[allen2Ind],sep='_'))
allen=cbind(allen1, allen2[,-1])
allen$ensembl=gtf$ensembl[match(allen$gene, gtf$gene)]
# read in hpa data. will end up not using this dataset as it does not correlate very well eith the existing dataset, but will use it because it does allow us to convert some more gene names to ENSG
normHpa=data.frame(fread('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/normal/hpa/rna_tissue_hpa.tsv'))
normHpa=normHpa[which(normHpa$Tissue=='cerebral cortex'),]
normHpa=normHpa[match(allen$gene, normHpa$Gene.name),]
allen$ensembl[which(is.na(allen$ensembl))]=normHpa$Gene[which(is.na(allen$ensembl))]
xli=allen$gene[which(is.na(allen$ensembl))]
out <- data.frame(queryMany(xli, scopes="symbol", fields="ensembl.gene", species="human", returnall=F, return.as='DataFrame'))
out$ensembl=as.character(lapply(out$ensembl, unlist))
out=out[-intersect(which(duplicated(out$query) | duplicated(out$query, fromLast=T)), grep('^[0-9]', out$X_id)),]
out=out[match(allen$gene, out$query),]
allen$ensembl[which(is.na(allen$ensembl))]=out$ensembl[which(is.na(allen$ensembl))]
allen=allen[,c(1,ncol(allen), seq(2,(ncol(allen)-1)))]
allen$ensembl[which(allen$ensembl=='NULL')]=NA
x=cor(allen[,-c(1,2)])
anno1=data.frame(as.factor(gsub('.*_', '', colnames(x))))
colnames(anno1)='Location'
rownames(anno1)=colnames(x)
pdf('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/plots/sample_cor_cnv_genes.pdf')
pheatmap(x,fontsize_col=5, fontsize_row=5, annotation_row=anno1, cutree_rows=7, cutree_cols=7)
dev.off()
xClust=cutree(hclust(dist(1-x)), 7)
allen=allen[,c(1,2, 18, 66,31)]
allen=allen[which(allen$ensembl %in% intersect(allen$ensembl, rnaOut$ENSG)),]
rnaOut=rnaOut[which(rnaOut$ENSG %in% intersect(allen$ensembl, rnaOut$ENSG)),]
allen=allen[-which(is.na(allen$ensembl)),]
allen=allen[match(rnaOut$ENSG, allen$ensembl),]
rnaOut$gene=allen$gene
setwd('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/SF10711_counts')
# write out SF10711 read counts
for(i in grep('SF10711', colnames(rnaOut))){
	name=paste0(colnames(rnaOut)[i], '.txt')
	fileOut=data.frame(ENSG=rnaOut$ENSG, Counts=rnaOut[,i])
	write.table(fileOut, file=name, sep='\t', row.names=F, quote=F, col.names=F)
}
setwd('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/normal_counts')
#  write out allen brain atlas read counts
for(i in grep('^S', colnames(allen))){
	name=paste0(colnames(allen)[i], '.txt')
	fileOut=data.frame(ENSG=allen$ensembl, Counts=allen[,i])
	write.table(fileOut, file=name, sep='\t', row.names=F, quote=F, col.names=F)
}
# write out correlation matrix
setwd('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/inputs/')
outCor$ENSG=rnaOut$ENSG[match(paste0(outCor$start, outCor$end, outCor$chr), paste0(rnaOut$start,rnaOut$end, rnaOut$chr))]
outCor=outCor[-which(is.na(outCor$ENSG)),]
outCor=outCor[which(outCor$ENSG %in% rnaOut$ENSG),]
outCorN=data.frame(ensembl=rnaOut$ENSG, hugo_gene=rnaOut$gene, kendall_t=rep(0, nrow(rnaOut)), pearson_r=rep(0, nrow(rnaOut)), spearman_r=rep(0, nrow(rnaOut)))
outCorN$kendall_t[which(outCorN$ensembl %in% outCor$ENSG)]=1
outCorN$pearson_r[which(outCorN$ensembl %in% outCor$ENSG)]=1
outCorN$spearman_r[which(outCorN$ensembl %in% outCor$ENSG)]=1
colnames(outCorN)[1]='Entrez_Gene_Id'
outCorN$Entrez_Gene_Id=gsub('ENSG0*', '', outCorN$Entrez_Gene_Id)
outCorN$hugo_gene=rnaOut$gene
write.table(outCorN, file='cor_all.tsv', sep='\t', row.names=F, quote=F)
#}}}
# RUN CNV-KIT WITH RNA INPUT
#{{{
cnvkit.py import-rna -f counts SF10711_counts/* -n normal_counts/* -c cor_all.tsv -g ensembl-gene-info.hg38.tsv --output ../intermediate/out-summary.tsv --output-dir ../intermediate 
cd ~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/intermediate
for ea in SF*cnr;	do
	name=$(sed 's/\..*/\.cns/g' <<< $ea)
	sem -j 15 cnvkit.py segment $ea -o ../outputs/$name 
done
#}}}
# Plot the data
#{{{
setwd('~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/outputs')
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
outPlot=melt(out[,-c(1,2,3,100)])
cnsEx=data.frame(fread("SF10711_9_1_1.cns"))
pdf('cnvs_all.pdf')
for(chr in seq(1,22)){
	tempPlot=outPlot[which(gsub(':.*', '', gsub('chr', '', outPlot$change))==chr),]
	tempCns=cnsEx[which(cnsEx$chromosome==chr),]
	print(ggplot(tempPlot, aes(x=variable, y=value, group=change)) + stat_smooth(aes(x=variable, y=value, color=change), method = "lm", formula = y~poly(x, 15), se = TRUE)+geom_point(aes(color=change))+theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4))+labs(x='Sections', y='logR', color='CNA', title='CNV-kit RNA infered CNA frequencies across sections', subtitle=paste0('Chromosome: ',chr, '; number of genes:', paste(tempCns$probes, collapse=','), ';Smoothing with 15-term polyn., with SE\n0.1 and -0.1 logR represent infered limit of detection from exome-sequencing'))+theme(legend.position = "bottom")+geom_hline(yintercept = .1)+geom_hline(yintercept = -.1))
}
dev.off()
outF=out[c(3,4,5,6,seq(15,20), 23,24, 31, 35,36, 38),-c(100,101)]
# to average, 15,16-7, 17,18-8, 19,20-9, include all called cnvs, even especially when splitting up arms of chromosome more unbiased
write.csv(outF, file='~/bdata/SF10711/cnv.analysis/rna.seq/cnvkit-rna/outputs/called_cnvs.csv')
#}}}
# Running Pyclone
#{{{
source('/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R')
library('stringr')
# reading amplicon sequencing data
amp=read.table('~/bdata/SF10711/amp.seq/csv_files/SF10711_variant_frequencies.csv', sep=',', header=T)
amp=amp[grep('SF10711_9_', amp$sample),]
names=paste('sf', str_pad(gsub('.*_', '', amp$sample), width=3, pad=0, side='left'),sep='')
genes_to_include=c('PTPRS', 'C10orf54', 'CARD14', 'GALNT2', 'THSD7B', 'LAMA2', 'MLX', 'PLEKHA5', 'MLL4', 'KRT17', 'OR8K3', 'INA', 'RYR2', 'CDKN1A', 'ATRX', 'TP53', 'TRIB2', 'TTN', 'TMCO4', 'PPIG', 'SSH1', 'PTPRZ1', 'ZSCAN10', 'HOXD4', 'OR11L1', 'IDH1', 'LRP2')
amp=amp[,c(1, grep(paste(genes_to_include, collapse='|'), colnames(amp)))]
amp=amp[grep('SF10711_9_', amp$sample),]
colnames(amp)=gsub('C10orf54', 'VSIR', colnames(amp)) 
colnames(amp)=gsub('C14orf37', 'ARMH4', colnames(amp))
colnames(amp)=gsub('MLL4', 'KMT2B', colnames(amp))
# which gene names to include
genes=unique(gsub('_.*', '', colnames(amp)[-1]))
genes[which(genes=='C10orf54')]='VSIR'
genes[which(genes=='C14orf37')]='ARMH4'
genes[which(genes=='MLL4')]='KMT2B'
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
# data.frame(cnvs_uN,meds=unlist(lapply(strsplit(cnvs_uN$Cellular.Fraction, ','),function(x) median(as.numeric(x))))) 
# cnvColDf is a dataframe that lists the CNVs across sections derived by FACETS
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
# excluding cnvs that do not align to the pattern preicted by the exome data
cnvKitOut=cnvKitOut[-4,]
cnvKitOut=data.frame(name=c('chr2pq_del', 'chr2q_del', 'chr11q_amp', 'chr17q_amp', 'chr18p_del'), cnvKitOut)
rownames(cnvKitOut)=cnvKitOut$name
cnvKitOut=cnvKitOut[,-seq(1,5)]
cnvKitOut=t(cnvKitOut)
cnvKitOut=cnvKitOut[which(rownames(cnvKitOut) %in% amp$sample),]
cnvKitOut=abs(cnvKitOut)
# read in GTF data so we can exact location of SNPS
gtf=fread('/home/shared/hg_align_db/GRCh37.38/gencode.v38lift37.annotation.gene_only.gtf', sep='\t')
gtf=gtf[,-c(2,3,6,8)]
colnames(gtf)=c('chr', 'start', 'end','strand', 'gene')
gtf$gene=gsub('";.*', '', gsub('.*gene_name "', '', gtf$gene))
gtf$gene[which(gtf$gene=='C10orf54')]='VSIR'
gtf$gene[which(gtf$gene=='C14orf37')]='ARMH4'
gtf$gene[which(gtf$gene=='MLL4')]='KMT2B'
gtf$mean=apply(gtf[,c(2,3)], 1, mean)
gtf=gtf[match(genes, gtf$gene),]
gtf$chr=gsub('chr', '', gtf$chr)
gtf=gtf[order(gtf$chr, gtf$start),]
gtf=cnAssign(cnv22, 's22')
gtf=cnAssign(cnv46, 's46')
gtf=cnAssign(cnv85, 's85')
gtf=cnAssign(cnv123,'s123')
gtf[which(gtf[,8]>3),8]=4
# major/minor copy number determined by closest exome sequencing section
sections=as.numeric(gsub('.*_', '', amp$sample))
sectionsDist=data.frame(abs(sections-22), abs(sections-46), abs(sections-85), abs(sections-123))
sectionsMin=apply(sectionsDist, 1, which.min)
setwd('~/bdata/SF10711/clonality/pyclone/inputs')
for(i in seq(1, nrow(amp))){
	mut=gsub('_Q.*', '', colnames(amp)[grep('depth', colnames(amp))])
	mut2=gsub('.*_', '', mut)
	mut=gsub('_.*', '', mut)
	mut[which(mut=='C10orf54')]='VSIR'
	mut[which(mut=='C14orf37')]='ARMH4'
	mut[which(mut=='MLL4')]='KMT2B'
	mut=paste(mut, mut2, sep='_')
	ref=amp[i,grep('\\.ref$', colnames(amp))]
	alt=amp[i,grep('\\.alt$', colnames(amp))]
	ind1=as.numeric(6+3*sectionsMin[i])
	ind2=as.numeric(5+3*sectionsMin[i])
	out=data.frame(mutation_id=mut, ref_counts=as.numeric(ref), var_counts=as.numeric(alt), normal_cn=rep(2, length(mut)), minor_cn=gtf[match(gsub('_.*', '', mut), gtf$gene),..ind1], major_cn=gtf[match(gsub('_.*', '', mut), gtf$gene),..ind2])
	colnames(out)=c('mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn')
	out$minor_cn[which(is.na(out$minor_cn))]=0
	out$major_cn[which(is.na(out$major_cn))]=2
	pseudosnps=data.frame(mutation_id=colnames(cnvKitOut), ref_counts=as.numeric(round(mean(out$ref_counts)*(1-cnvKitOut[i,]))), var_counts=as.numeric(round(mean(out$ref_counts)*(cnvKitOut[i,]))), normal_cn=rep(2, ncol(cnvKitOut)), minor_cn=rep(0, ncol(cnvKitOut)), major_cn=rep(1, ncol(cnvKitOut)))
	out=rbind(out, pseudosnps)
	write.table(out, file=paste(names[i], 'tsv', sep='.'), sep='\t', row.names=F)
}
out_com=paste('~/bdata/SF10711/clonality/pyclone/inputs/', names, '.tsv', sep='',collapse=',')
out_com=paste('PyClone run_analysis_pipeline --in_files {', out_com, '} ','--tumour_contents {',paste(amp$TP53_G338T_freq, collapse=','), '} ', '--working_dir ~/bdata/SF10711/clonality/pyclone/outputs', sep='')
setwd('/mnt/bdata/patrick/SF10711/clonality/pyclone/')
write.table(out_com, file='pyclone_command.txt', quote=F, row.names=F, col.names=F)
#}}}
# MAKE SURE YOU ARE NOT USING X-FORWARDING and run in SCREEN or TMUX...
# run PyClone using printed output command
{{{
cd /mnt/bdata/patrick/SF10711/clonality/pyclone
conda activate python2
chmod +x pyclone_command.txt
./pyclone_command.txt
conda deactivate
}}}
# convert output to cITUP format
# test for clonal copy number alterations that were for corrected using Pyclone
#{{{
setwd('~/bdata/SF10711/clonality/pyclone/outputs/tables')
clone=read.csv('loci.tsv', sep='\t')
clust=read.csv('cluster.tsv', sep='\t')
clone=clone[,c(1,3)]
cloneOut=unique(apply(clone,1, paste, collapse='-'))
cloneClone=gsub('.*-', '', cloneOut)
cloneName=gsub('_.*', '', cloneOut)
setwd('~/bdata/SF10711/clonality/pyclone/outputs/trace')
i=1
for(ea in list.files()[grep('.cellular_prevalence.tsv.bz2', list.files())]){
inFile=read.table(ea, sep='\t', header=T)[-seq(1,100),]
inFile=c(apply(inFile, 2, mean))
if(i==1){outFile=data.frame(matrix(nrow=length(inFile), ncol=length(grep('.cellular_prevalence.tsv.bz2', list.files()))))}
outFile[,i]=inFile
colnames(outFile)[i]=gsub('\\..*', '', ea)
i=i+1
}
rownames(outFile)=names(inFile)
amp=t(outFile)
colnames(amp)=gsub('C10orf54', 'VSIR', colnames(amp))
colnames(amp)=gsub('C14orf37', 'ARMH4', colnames(amp))
colnames(amp)=gsub('MLL4', 'KMT2B', colnames(amp))
colnames(amp)=gsub('_.*', '', colnames(amp))
adjAmp=amp
amp=read.table('~/cluster/big-data2/Sam/Loglio_project/SF10711_astrocytoma/AmpSeq_data/AmpSeq_analysis/Aligned_files_bwa_parallelized/Mutation_variant_frequency_data/SF10711_variant_frequencies.csv', sep=',', header=T)
amp=amp[grep('SF10711_9_', amp$sample), grep('freq', colnames(amp))]
colnames(amp)=gsub('_.*', '', colnames(amp))
colnames(amp)=gsub('C10orf54', 'VSIR', colnames(amp))
colnames(amp)=gsub('C14orf37', 'ARMH4', colnames(amp))
colnames(amp)=gsub('MLL4', 'KMT2B', colnames(amp))
colnames(amp)=gsub('_.*', '', colnames(amp))
adjAmp=adjAmp[,match(colnames(amp), colnames(adjAmp))]
diff=adjAmp[,]-amp[,]
out=apply(diff, 2, mean)
ampM=apply(amp, 2, mean)
adjAmpM=apply(adjAmp,2, mean)
gtf$diff=out[match(gtf$gene, names(out))]
gtf$adjamp=adjAmpM[match(gtf$gene, names(adjAmpM))]
gtf$amp=ampM[match(gtf$gene, names(ampM))]
gtf$difScale=gtf$adjamp/gtf$amp
# convert output to cITUP format
setwd('~/bdata/SF10711/clonality/pyclone/outputs/tables')
clone=read.csv('loci.tsv', sep='\t')
clone=clone[,c(1,3)]
cloneOut=unique(apply(clone,1, paste, collapse='-'))
cloneClone=gsub('.*-', '', cloneOut)
cloneName=gsub('-.*', '', cloneOut)
setwd('~/bdata/SF10711/clonality/pyclone/outputs/trace')
i=1
for(ea in list.files()[grep('.cellular_prevalence.tsv.bz2', list.files())]){
	inFile=read.table(ea, sep='\t', header=T)[-seq(1,100),]
	inFile=c(apply(inFile, 2, mean))
	if(i==1){outFile=data.frame(matrix(nrow=length(inFile), ncol=length(grep('.cellular_prevalence.tsv.bz2', list.files()))))}
	outFile[,i]=inFile
	colnames(outFile)[i]=gsub('\\..*', '', ea)
	i=i+1
}
rownames(outFile)=names(inFile)
amp=t(outFile)
colnames(amp)=gsub('C10orf54', 'VSIR', colnames(amp))
colnames(amp)=gsub('C14orf37', 'ARMH4', colnames(amp))
colnames(amp)=gsub('MLL4', 'KMT2B', colnames(amp))
indChange=match(colnames(amp), cloneName)
cloneName=cloneName[indChange]
cloneClone=cloneClone[indChange]
amp[(amp<0)]=1E-4
cloneClone=as.numeric(cloneClone)-1
setwd('~/bdata/SF10711/clonality/pyclone_cITUP/inputs')
write.table(t(amp), file='freq.txt', quote=F, row.names=F, col.names=F)
write.table(cloneClone, file='clusters.txt', quote=F, row.names=F, col.names=F)
#}}}
# in bash shell run cITUP data
#{{{
cd ~/bdata/SF10711/clonality/pyclone_cITUP
conda activate pyclone
run_citup_qip.py --maxjobs 10 --min_nodes 5 --max_nodes 5 --nocleanup ./inputs/freq.txt ./inputs/clusters.txt ./outputs/results.h5
conda deactivate
#}}}
# Plot tree from cITUP output
#{{{
source('/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R')
citupToTree(citDir='~/bdata/SF10711/clonality/pyclone_cITUP/outputs', snpDir='~/bdata/SF10711/clonality/pyclone/outputs/tables', outDir='~/bdata/SF10711/clonality/pyclone_cITUP/trees')
#}}}
