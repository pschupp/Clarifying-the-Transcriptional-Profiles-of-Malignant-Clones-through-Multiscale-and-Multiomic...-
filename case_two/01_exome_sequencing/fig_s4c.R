# dependencies and custom function
# {{{
library('circlize')
library('ComplexHeatmap')
library('stringr')
library('svglite')
library('extrafont')
library('gridBase')
library('gridtext')
vec2col=function(col1, col2, n, vec, col3=NA){
	if(class(vec)=='numeric'){
		naNumb=length(which(is.na(vec)))
		if(is.na(col3)){mypal <- colorRampPalette(c(col2,col1))(n)} else{mypal <- colorRampPalette(c(col2,col3,col1))(n)}
		vec[which(vec==Inf)]=min(vec[-which(vec==Inf)])
		if(max(vec,na.rm=T)>1.1){names(mypal)=seq(0, 100, length=n)}else{names(mypal)=seq(0, 1, length=n)}
		out=rep('#5b5b5b', length(vec))
		i=1
		for(ea in vec){
			if(is.na(ea)){i=i+1;next}
			ind=which.min(abs(as.numeric(names(mypal))-ea))
			out[i]=mypal[ind]
			i=i+1
		}
	}
	if(class(vec)=='factor'){
		out=rep(col1, length(vec))
		out[which(vec=='No')]=col2
		out[which(vec==0|is.na(vec)|vec==F|vec==FALSE)]='#5b5b5b'
	}
	out
}
# }}}
# read in exome, amplicon, and TCGA mutation data
# {{{
exome=read.table('~/@patrick/SF10711/figures/tables/Summary_of_nonsynonamous_exonic_mutations_SF10711.txt', header=T, sep='\t')
amp=read.table('Piece_9_High_confidence_mutations_>200reads>10xfreq_in_blood.csv', sep=',', header=T)
amp=amp[,c(1, grep('*_freq', colnames(amp)))]
ampfull=read.csv('SF10711_variant_frequencies.csv')
ampfull=ampfull[,c(1, grep('*_freq', colnames(ampfull)))]
rna=data.frame(fread('ERCC_K10_RUVg_normalized_cpm_values_for_RNAseq.csv', sep=',', header=T))
rna=rna[,-c(2,3,4,5,6)]
setwd('..')
change=gsub('ins.*','ins',as.character(unlist(sapply(exome$protein, function(x) strsplit(x, ',')[[1]][1]))))
change[which(is.na(change))]=exome$nucleotide[which(is.na(exome$nucleotide))]
change[which(is.na(change))]='intronic'
dat=data.frame(Label=paste(exome$gene, change, sep='_'))
ampDetec=rep(NA, length(dat$Label))
ampDetec[which(gsub('_.*','', dat$Label) %in% gsub('_.*', '', colnames(ampfull)[-1]))]='No'
ampDetec[which(gsub('_.*','', dat$Label) %in% gsub('_.*', '', colnames(amp)[-1]))]='Yes'
mut=read.csv("~/cluster/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/TCGA_LGG_data_030615/LGG_FINAL_ANALYSIS.aggregated.capture.tcga.uuid.curated.somatic.maf.csv")
bio1=read.csv('~/cluster/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/TCGA_LGG_data_030615/Clinical/Biotab/nationwidechildrens.org_biospecimen_cqcf_lgg.txt', sep='\t') 
mut$TumorSample=unlist(lapply(strsplit(mut$Tumor_Sample_Barcode,'-'), function(x) paste(x[1:3], collapse='-')))
bio1=bio1[which(bio1$histological_type=='Astrocytoma'),]
mut=mut[which(mut$TumorSample %in% bio1$bcr_patient_barcode),]
mutOut=data.frame(hugoGene=mut$Hugo_Symbol, type=mut$Variant_Classification, chr=mut$Chromosome, start=mut$Start_position, end=mut$End_position, sample=mut$TumorSample, alt=mut$t_alt_count, ref=mut$t_ref_count)
mutOut=data.frame(mutOut, vaf=mutOut$alt/(mutOut$ref+mutOut$alt))
mutTot=aggregate(data.frame(mutOut$vaf), by=list(mutOut$hugoGene), FUN="mean")
colnames(mutTot)=c('Gene', 'vaf')
# }}}
# tidy some labels
# {{{
dat=data.frame(exome[,c(1,2,3,4,5)], section22=exome$Recurrence1_9.1.22.var.freq, section46=exome$Recurrence1_9.1.46.var.freq, section85=exome$Recurrence1_9.2.85.var.freq, section123=exome$Recurrence1_9.2.123.var.freq, ampseq=ampDetec, sanger=(gsub('_.*', '', dat$Label) %in% c('IDH1', 'PKD1L1', 'MLH1', 'CASP1', 'ATRX', 'NXRN1', 'AXDND1')), TCGA=mutTot$vaf[match(gsub('_.*', '', dat$Label), mutTot$Gene)])
rownames(rna)=rna$Geneid
rna=rna[,-1]
rnaMean=rowMeans(rna)
rnaPerc=(order(rnaMean,decreasing=T)/length(rnaMean))*100
names(rnaPerc)=names(rnaMean)
dat=data.frame(dat, MeanExprPerc=rnaPerc[match(dat$gene, names(rnaPerc))])

setwd('~/@patrick/SF10711/figures/fig3/wheelplot/SF10711/inputs')
exome=read.table('Summary_of_nonsynonamous_exonic_mutations_SF10711.txt', header=T, sep='\t')
out=data.frame(chr=gsub('chr','', exome$chr), pos=exome$position,  change=apply(exome[,c(4,5)], 1, paste, collapse='/'), strand=rep('+', nrow(exome)))
allA=exome[,4]
allA[allA=='-']=''
out=data.frame(out, end=out$pos+(nchar(allA)-1))
out=data.frame(out$chr, start=out$pos, end=out$end, out$change, out$strand)
out$end[grep('^-/', out$out.change)]=out$start[grep('^-/', out$out.change)]-1
write.table(out, col.names=F, row.names=F, sep=' ', file='SF10711.positions.for.vep.txt', quote=F)
# }}}
# use vep to convert to h38 and get official hgvs names
# {{{
# bash
vep	-i SF10711.positions.for.vep.txt \
    -o SF10711.vep.output.new \
    --total_length \
    --cache --dir_cache /home/shared/vep_cache/ \
    --fasta /home/shared/vep_cache/fasta/GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa \
    --hgvs --hgvsg --shift_hgvs 1 --symbol --refseq --fork 10 --pick --offline  --numbers  --assembly GRCh37 \
    --tab  --fields  "Uploaded_variation,SYMBOL,Location,STRAND,Allele,Gene,Feature,Consequence,SIFT,POLYPHEN,IMPACT,HGVSc,HGVSp,HGVSg,EXON,INTRON,cDNA_position,CDS_position,Protein_position" \
    --force
setwd('~/@patrick/SF10711/figures/liftover_GRCh37_to_38')
vep=fread('SF10711.vep.output', header=T, sep='\t', skip='#Uploaded_variation')
dat=data.frame(dat$gene, dat$chr, dat$position, vep$STRAND, dat$ref_allele, dat$alt_allele, vep$Feature, vep$Consequence, vep$HGVSc, vep$HGVSp, dat$section22, dat$section46, dat$section85, dat$section123, dat$sanger, dat$ampseq, dat$MeanExprPerc, dat$TCGA)
colnames(dat)=c('Gene', 'Chromosome', 'Position', 'Strand', 'Reference_allele', 'Alternate_allele', 'Refseq_transcript_ID', 'Consequence', 'HGVSc', 'HGVSp', 'Section_22_VAF', 'Section_46_VAF', 'Section_85_VAF', 'Section_123_VAF', 'Verified_sanger', 'Verified_ampseq', 'Mean_expression_percentile', 'TCGA_astrocytoma_vaf')
dat=dat[order(apply(dat[,grep('Section.*VAF$', colnames(dat))], 1, mean, na.rm=T), decreasing=T),]
write.table(dat, file='SF10711_wheelplot_table.tsv', sep='\t', quote=F, row.names=F)
setwd('~/@patrick/SF10711/figures/liftover_GRCh37_to_38')
input=fread('SF10711.positions.for.vep.txt')
output=fread('SF10711.vep.output')
exome=read.table('SF10711_wheelplot_table.tsv', header=T, sep='\t')
write.table(data.frame(paste('chr',input$V1, sep=''), input[,2], input[,3]+1, input[,5], output$SYMBOL), sep='\t', col.names=F, row.names=F, quote=F, file='positions_for_liftover.tsv')
PicardCommandLine CreateSequenceDictionary  R=/home/shared/hg_align_db/GRCh37.38/GRCh37.primary_assembly.genome.fa \
                                            O=/home/shared/hg_align_db/GRCh37.38/GRCh37.primary_assembly.genome.seqDict.Picard.dict \
                                            GENOME_ASSEMBLY=GRCh37
cat /home/shared/hg_align_db/GRCh37.38/GRCh37.primary_assembly.genome.seqDict.Picard.dict positions_for_liftover.tsv  > test.txt
mv test.txt positions_for_liftover.tsv
# Interval list format per (https://gatk.broadinstitute.org/hc/en-us/articles/360040096472-IntervalListTools-Picard-), 
# tab seperated [chromsome] [start] [stop]  [strand]    [name]
PicardCommandLine LiftOverIntervalList  INPUT=positions_for_liftover.tsv \
                                        OUTPUT=GRCh38_positions.tsv \
                                        SEQUENCE_DICTIONARY=/home/shared/hg_align_db/GRCh37.38/GRCh37.primary_assembly.genome.seqDict.Picard.dict \
                                        REJECT=reject.txt \
                                        CHAIN=/home/shared/hg_align_db/liftover/hg19ToHg38.over.primary.chain \
                                        MIN_LIFTOVER_PCT=0.95

exome=read.table('~/@patrick/SF10711/figures/tables/sf10711_wheelplot_table.tsv', header=T, sep='\t', fill=T)
exome=exome[order(exome$Chromosome, exome$Position),]
output=fread('~/@patrick/SF10711/figures/liftover_GRCh37_to_38/SF10711.vep.output')
output$chr=gsub(':.*', '', output$Location)
output$pos=as.numeric(gsub('-.*', '', gsub('.*:', '', output$Location)))
output=output[order(output$chr, output$pos),]
lift=fread('~/@patrick/SF10711/figures/liftover_GRCh37_to_38/GRCh38_positions.tsv', sep='\t', skip=85)
data.frame(exome$Gene, output$SYMBOL, lift$V5)
exome=exome[order(exome$Chromosome, exome$Position),]
convert=data.frame(old=exome$Gene, new=paste(output$SYMBOL, gsub('.*:p\\.', '', output$HGVSp), sep=' '))
convert[grep('-', convert$new),2]=paste(output$SYMBOL[grep('-', convert$new)], gsub('_', ' ', gsub(',.*', '', output$Consequence[grep('-', convert$new)])))
convert$new=gsub('splice.*', 'splice', convert$new)
convert$new=gsub('Leu', 'L', convert$new)
convert$new=gsub('Lys', 'K', convert$new)
convert$new=gsub('Ile', 'I', convert$new)
convert$new=gsub('Pro', 'P', convert$new)
convert$new=gsub('Phe', 'F', convert$new)
convert$new=gsub('Arg', 'R', convert$new)
convert$new=gsub('His', 'H', convert$new)
convert$new=gsub('Ser', 'S', convert$new)
convert$new=gsub('Cys', 'C', convert$new)
convert$new=gsub('Val', 'V', convert$new)
convert$new=gsub('Thr', 'T', convert$new)
convert$new=gsub('Tyr', 'Y', convert$new)
convert$new=gsub('Met', 'M', convert$new)
convert$new=gsub('Trp', 'W', convert$new)
convert$new=gsub('Ala', 'A', convert$new)
convert$new=gsub('Asp', 'D', convert$new)
convert$new=gsub('Asn', 'N', convert$new)
convert$new=gsub('Gly', 'G', convert$new)
convert$new=gsub('Gln', 'N', convert$new)
convert$new=gsub('Glu', 'E', convert$new)
write.table(convert, row.names=F, quote=F, sep=',', file='sf10711_muts_liftover.csv')

exome$Position=lift$V2
exome$Gene=output$SYMBOL
exome$Consequence=output$Consequence
exome$HGVSc=output$HGVSc
exome$HGVSp=output$HGVSp
write.table(exome, row.names=F, quote=F, sep='\t', file='~/@patrick/SF10711/figures/tables/sf10711_wheelplot_table.tsv')
# }}} 
# tidy labels
# {{{
setwd('~/@patrick/SF10711/figures')
dat=read.table('~/@patrick/SF10711/figures/tables/sf10711_wheelplot_table.tsv', sep='\t', header=T, row.names=NULL)
dat=dat[which(dat$Verified_ampseq=='Yes'),]
geneSelec=c('PTPRS', 'VSIR', 'CARD14', 'GALNT2', 'THSD7B', 'LAMA2', 'MLX', 'PLEKHA5', 'KMT2B', 'KRT17', 'OR8K3', 'INA', 'RYR2', 'CDKN1A', 'ATRX', 'TP53', 'TRIB2', 'TTN', 'TMCO4', 'PPIG', 'SSH1', 'PTPRZ1', 'ZSCAN10', 'HOXD4', 'IDH1', 'LRP2', 'RUFY1')
dat=dat[-ncol(dat),]
dat$Mutation.distribution=apply(dat[,seq(11,14)], 1, function(x) paste(c(22,46,85,123)[-which(is.na(x))], collapse=','))
dat$Mutation.distribution[apply(dat[,seq(11,14)], 1,function(x) length(which(!(is.na(x))))==4)]='22,46,85,123'
write.table(dat,file='~/@patrick/SF10711/figures/tables/sf10711_wheelplot_table.tsv', sep='\t', quote=F, row.names=F)
dat=dat[which(dat$Gene %in% geneSelec),]
write.table(dat, row.names=F, quote=F, file='~/@patrick/SF10711/figures/tables/sf10711_wheelplot_table_publication_snps.tsv', sep='\t')
# }}}
# create wheelplot 
# {{{
dat = fread('~/@patrick/SF10711/figures/tables/sf10711_wheelplot_table_publication_snps.tsv', sep='\t', data.table = FALSE)
dat=data.frame(dat, Label=paste(dat$Gene, gsub('.*:', '', dat$HGVSp), sep=':'))
dat$Label[which(dat$HGVSp=='-')]=paste(dat$Gene[which(dat$HGVSp=='-')], gsub(',.*', '', dat$Consequence[which(dat$HGVSp=='-')]), sep=':')
dat$Label=gsub('_variant', '', dat$Label)
dat=dat[order(apply(dat[,grep('^Section', colnames(dat))], 1, mean, na.rm=T),decreasing=T),]
nMax=max(unlist(lapply(dat$Label, nchar)))
new=floor(ceiling((nrow(dat))/.75)/2)
dat$Mean_expression_percentile[is.na(dat$Mean_expression_percentile)]=0
oldRow=35
naDat=data.frame(matrix(NA, nrow=ceiling(ceiling(new)/2), ncol=ncol(dat)))
colnames(naDat)=colnames(dat)
naDat$Label=make.names(naDat$Label, unique=T)
dat$Section_22_VAF[is.na(dat$Section_22_VAF)]=0
dat$Section_46_VAF[is.na(dat$Section_46_VAF)]=0
dat$Section_85_VAF[is.na(dat$Section_85_VAF)]=0
dat$Section_123_VAF[is.na(dat$Section_123_VAF)]=0
dat=data.frame(rbind(dat, naDat))
dat[,c(1,2,4,5,6,7,8,9,10,15,16,19)]=lapply(dat[,c(1,2,4,5,6,7,8,9,10,15,16,19)], as.factor)
dat$Label=factor(dat$Label, levels=dat$Label)
nameCol=rep('#5b5b5b', length(dat$Label))
nameCol[which(dat$ampseq=='Yes')]='black'
names(nameCol)=dat$Label
dat$Label=gsub(':(p\\.)?', ' ', dat$Label)
dat$Label=gsub('Leu', 'L', dat$Label)
dat$Label=gsub('Lys', 'K', dat$Label)
dat$Label=gsub('Ile', 'I', dat$Label)
dat$Label=gsub('Pro', 'P', dat$Label)
dat$Label=gsub('Phe', 'F', dat$Label)
dat$Label=gsub('Arg', 'R', dat$Label)
dat$Label=gsub('His', 'H', dat$Label)
dat$Label=gsub('Ser', 'S', dat$Label)
dat$Label=gsub('Cys', 'C', dat$Label)
dat$Label=gsub('Val', 'V', dat$Label)
dat$Label=gsub('Thr', 'T', dat$Label)
dat$Label=gsub('Tyr', 'Y', dat$Label)
dat$Label=gsub('Met', 'M', dat$Label)
dat$Label=gsub('Trp', 'W', dat$Label)
dat$Label=gsub('Ala', 'A', dat$Label)
dat$Label=gsub('Asp', 'D', dat$Label)
dat$Label=gsub('Asn', 'N', dat$Label)
dat$Label=gsub('Gly', 'G', dat$Label)
dat$Label=gsub('Gln', 'N', dat$Label)
dat$Label=gsub('Glu', 'E', dat$Label)
dat$Label=gsub('splice_.*', 'splice', dat$Label)
dat$Label=factor(dat$Label, levels=dat$Label)
svg(paste0(WD, paste('SF10711.wheelplot', oldRow, 'svg', sep='.')), height=10, width=7.2)
plot.new()
circle_size = unit(1, "snpc") 
pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
par(omi = gridOMI(), new = TRUE)
circos.clear()
circos.par(points.overflow.warning=F, track.height=0.105, start.degree=90, circle.margin=1, track.margin=c(0.006, 0.006))
circos.initialize(as.factor(dat$Label), xlim=c(-1,1))
circos.track(dat$Label,bg.col=vec2col('#4daf4a', 'white', 100, dat$Section_22_VAF), ylim=c(-1,1), bg.lwd=1.5, panel.fun = function(x, y) {
	xlim = get.cell.meta.data("xlim")
	ylim = get.cell.meta.data("ylim")
	colVec=nameCol[which(names(nameCol)==CELL_META$sector.index)]
	labelVec=CELL_META$sector.index
	circos.text(x=CELL_META$xcenter,y=ylim[1]+mm_y(6), col=colVec ,labels=,labelVec, facing = "clockwise", niceFacing = TRUE, family="NimbusSan", font=2, cex=1.2, adj=c(0,0.5))
})
circos.track(dat$Label,bg.col=vec2col('#4daf4a', 'white', 100, dat$Section_46_VAF), ylim=c(-1,1), bg.lwd=1.5)
circos.track(dat$Label,bg.col=vec2col('#4daf4a', 'white', 100, dat$Section_85_VAF), ylim=c(-1,1), bg.lwd=1.5)
circos.track(dat$Label,bg.col=vec2col('#4daf4a', 'white', 100, dat$Section_123_VAF), ylim=c(-1,1), bg.lwd=1.5)
colTemp=vec2col('black', 'white', 100, dat$Verified_ampseq)
circos.track(dat$Label,bg.col=colTemp, ylim=c(-1,1), bg.lwd=1.5)
circos.track(dat$Label,bg.col=vec2col('#377eb8', 'white', 100, dat$TCGA_astrocytoma_vaf), ylim=c(-1,1),force.ylim=T, bg.lwd=1.5)
circos.track(dat$Label,bg.col=vec2col('#e31a1c', 'white', 100, dat$Mean_expression_percentile), ylim=c(-1,1), bg.lwd=1.5)
rect(-2, -.0,-.015,2, col='white', border=NA)
textAdj=0.13
textStart=.98
text(x=-0.7, y=textStart-0*textAdj, labels=str_pad('22', 32, 'left'),family='NimbusSan', cex=1.1, col='#4daf4a')
text(x=-1.05, y=textStart-1.4*textAdj, labels=str_pad('Exome VAF      ', 32, 'left'),family='NimbusSan', cex=1.1, col='#4daf4a')
text(x=-3,y=textStart-1.4*textAdj, labels=str_pad('{', 32, 'left'),family='NimbusSan', cex=5, col='#4daf4a')
text(x=-0.7, y=textStart-1*textAdj, labels=str_pad('46', 32, 'left'),family='NimbusSan', cex=1.1, col='#4daf4a')
text(x=-.7, y=textStart-2*textAdj, labels=str_pad('86', 32, 'left'),family='NimbusSan', cex=1.1, col='#4daf4a')
text(x=-0.7, y=textStart-3*textAdj, labels=str_pad('123', 32, 'left'),family='NimbusSan', cex=1.1, col='#4daf4a')
text(x=-0.85, y=textStart-4*textAdj, labels=str_pad('Amp-seq verified', 32, 'left'),family='NimbusSan', cex=1.1, col='black')
text(x=-1, y=textStart-5*textAdj, labels=str_pad('Gene-based VAF, TCGA astro.', 32, 'left'),family='NimbusSan', cex=1.1, col='#377eb8')
text(x=-0.98, y=textStart-6*textAdj, labels=str_pad('Mean expression percentile', 32, 'left'),family='NimbusSan', cex=1.1, col='#e31a1c')
legendTCGA=Legend(at=c(0,1), col_fun=colorRamp2(c(1,0), c('#377eb8', 'white')), title='Gene-based VAF\nTCGA astrocytomas', border=F, title_position='topcenter', title_gp=gpar(fontsize=13, fontface="bold", fontfamily='NimbusSan'), labels_gp=gpar(fontsize=13, fontface="bold", fontfamily='NimbusSan'))
legendExome=Legend(at=c(0,1),col_fun=colorRamp2(c(1,0), c('#4daf4a', 'white')), title='Exome\nVAF', border=F, title_position='topcenter', title_gp=gpar(fontsize=13, fontface="bold", fontfamily='NimbusSan'), labels_gp=gpar(fontsize=13, fontface="bold", fontfamily='NimbusSan'))
legendExpr=Legend(at=c(0,100),col_fun=colorRamp2(c(100,0), c('#e31a1c','white')), title='Mean expr.\npercentile', border=F, title_position='topcenter', title_gp=gpar(fontsize=13, fontface="bold", fontfamily='NimbusSan'), labels_gp=gpar(fontsize=13, fontface="bold", fontfamily='NimbusSan'))
upViewport()
draw(packLegend(legendExome,legendTCGA,legendExpr, direction="horizontal", gap=unit(50, 'mm')), y=unit(1,"npc")-circle_size, just="top")
dev.off()
# }}}
