## Unzip the data from /SJS_01 and copy the fastq file to expanded_fastq_first_lane:
# {{{
for f in *.gz; do
STEM=$(basename "${f}" .gz)
gunzip -c "${f}" > /big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/expanded_fastq_first_lane/"${STEM}"
done
# }}}
## Now I would like to run Fastqc on the files to see how many reads I am getting for each section and what the quality of the reads looks like
# {{{
## Add fastqc to my PATH
## Run fastqc on one of the files prior to adaptor trimming
/fast-data/bin/jdk1.8.0_31/bin/java fastqc SF10711_9-1-10_S22_L004_R1_001.fastq --outdir ../FastQC_output
for file in *.fastq
do
fastqc "$file" --outdir ../FastQC_output
done
## Now I would like to extact the coverage for each of the samples from the summary statistics file. To do this I first need to unzip each of the zip files.
for file in *.zip
do
done
## Now I will read in the fastqc_data.txt file into R and extract the relevant information and put this into a table.
R
setwd("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/FastQC_output")
list1=list.dirs(getwd())
## Subset the list down to just folders with the fastqc files in.
list1=list1[-(grep("/Images",list1,fixed=TRUE))]
list1=list1[-(grep("/Icons",list1,fixed=TRUE))]
list1=list1[-(grep("Undetermined",list1,fixed=TRUE))]
list1=list1[-1]
length(list1)
#[1] 96
## I now have 96 folders which I need to read the information in and extract the statistics from.
mat1=matrix(nrow=96,ncol=6)
colnames(mat1)=c("sample_name","total_reads","Pcent_GC","Per_base_sequence_quality","Overrepresented_sequences","Pcent_of_sequence")
for (i in 1:96){
    read1=read.delim(paste(list1[i],"fastqc_data.txt",sep="/"))
    mat1[i,1]<- sapply(strsplit(as.character(read1[3,2]),split="_S",fixed=TRUE),"[",1)
    mat1[i,2]<-as.character(read1[6,2])
    mat1[i,3]<-as.character(read1[9,2])
    mat1[i,4]<-as.character(read1[11,2])
    
    if(length(grep("TruSeq Adapter",read1[,2]))==1){
    mat1[i,5]<-as.character(read1[grep("TruSeq Adapter",read1[,2]),2])
    mat1[i,6]<-as.character(read1[grep("TruSeq Adapter",read1[,2]),1])
    }else{
        mat1[i,5]<-paste(as.character(read1[grep("TruSeq Adapter",read1[,2]),2]), sep="_", collapse=":")
        mat1[i,6]<-paste(as.character(read1[grep("TruSeq Adapter",read1[,2]),1]), sep="_", collapse=":")
    }
}
mat1=as.data.frame(mat1)
section.no=sapply(strsplit(as.character(mat1[,1]),split="-",fixed=TRUE),"[",3)
section.no=as.numeric(section.no)
mat1=data.frame(sample_name=mat1[,1],section.no=section.no,mat1[,2:dim(mat1)[2]])
mat1=mat1[order(mat1$section.no),]
mat1$total_reads=as.character(mat1$total_reads)
mat1$total_reads=as.numeric(mat1$total_reads)
## What is the range of reads
range(mat1$total_reads)
#[1] 3027871 5321752
## If I were to get the same number of reads for the 7 additional lanes then this would give me coverage of between 24.2 and 42.6 million reads per sample:
8*3027871
#[1] 24222968
8*5321752
#[1] 42574016
## The difference between the highest coverage sample and the lowest coverage sample is 1.76 fold which is not a huge difference.
5321752/3027871
#[1] 1.757589
setwd("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/FastQC_output/Analysis_in_R")
pdf("histogram_of_reads.pdf")
hist(mat1$total_reads)
dev.off()
## The mean coverage is around 4 million reads per sample.
mean(mat1$total_reads)
#[1] 4041149
## Write this table out to disk
write.table(mat1,file="FastQC_summary_for_each_sample_SJS-01_first_lane.csv",col.names=TRUE,row.names=FALSE,sep=",")
# }}}
# trim galore
# {{{
# using a q20 cutoff, running fastqc and having a minimum length of 20 to retain the read.
trim_galore -q 20 --fastqc --length 20 -o ./trimmed --paired ${READ1[$i]} ${READ2[$i]}
# }}}
# generate human genome + ERCC
# {{{
cat $hg19list /big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/ERCC92_reference_genome/ERCC92.fa > ../combined_hg19_ERCC92_reference/combined_hg19_ERCC92.fa
## Now combine the gtf files.
cat /big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/hg19_reference_genome/annotation_files/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf /big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/ERCC92_reference_genome/ERCC92.gtf > ../combined_hg19_ERCC92_reference/combined_hg19_ERCC92.gtf
## Now build a bowtie2 index from the new combined fasta
bowtie2-build combined_hg19_ERCC92.fa ./Bowtie2Index/combined_hg19_ERCC92
# }}}
# alignment with tophat
# {{{
for file in *.fq
do
    tophat -p 32 -G combined_hg19_ERCC92.gtf --library-type=fr-firststrand -o ../stranded_hg19_ERCC_combined_Tophat_output_trimmed_reads_bowtie2/"${file/_S*_L*_R1_001_trimmed.fq/}"_thout combined_hg19_ERCC92 "$file"
done
}}}
# feature counts
# {{{
## The above code gives me a count on an individual gene basis but I am not sure that this makes sense. As far as I understand this it is counting all reads from all exons within a group. The group in this case is gene_id so it is counting the reads from all exons from a particular gene. Biologically this doesn't make sense because the predominant transcript doesn't necessarily have all exons so the overall result doesn't really make sense. What would probably make more sense is to group by transcript ID and then pick out the major transcript by its ID.
featureCounts -T 5 -a /big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/ERCC92_reference_genome/ERCC92.gtf -t exon -g gene_id -o /big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/featurecounts_output/Trimmed_reads_ERCC_controls_combined_count_by_gene_id.txt $filelist3
# }}}
# ERCC analysis
# {{{
reads=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/SF10711_mapped_to_hg19_and_ERCC92_stranded_round2_SampleNetworks/All_10-56-41/SF10711_mapped_to_hg19_and_ERCC92_stranded_round2_All_89_outliers_removed.csv")
dim(reads)
#[1] 22933    95
## Get the RNAseq samples in sequential order
readsSub=reads[,7:95]
toOrder=sapply(strsplit(as.character(colnames(readsSub)),split="_",fixed=TRUE),"[",4)
toOrder=as.numeric(toOrder)
readsSub=readsSub[,order(toOrder)]
reads=data.frame(reads[,1:6],readsSub)
## I would like to calculate the RLE this is the relative log expression. log-ratio of a read count to the median count across samples.
Genemedian=apply(reads[,7:95], 1, median)
RLE=apply(reads[,7:95],2,function(x) log2(x/Genemedian))
rownames(RLE)=reads$Geneid
## Looking at the paper by Davide Risso which looked at the ERCC controls they made box plots of un-normalized relative log2 expression (RLE) to see how even each of the libraries are.
RLE1=lapply(c(1:89), function(x){RLE[,x]})
names(RLE1)=colnames(RLE)
pdf(file="boxplot of normalized libraries RLE hg19 ERCC92 trimmed stranded.pdf",width=9,height=5)
par(mar=c(5,5,4,2))
boxplot(RLE1,notch=T,axes=FALSE,ylab="RLE",xlab="library",main="Relative log expression across serial sections of SF10711",cex=0.3,lwd=0.3)
# yaxis
axis(2, seq(-6,4, 2))
# xaxis
axis(1,1:length(RLE1),gsub("SF10711_","",names(RLE1),fixed=TRUE),las=2,cex.axis=0.5)
dev.off()
rm(list=ls())
## Now try normalization using RUVg and the ERCC controls as standards. This was desribed by Davide Risso and appears to be the best method for RNAseq normalization.
library(RUVSeq)
## Read in the unnormalized data
dat1=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/RNAseq_non-normalized-hg19-ERCC92_trimmed_reads_stranded.csv")
## Now create a matrix from the gene expression data
datExpr=as.matrix(dat1[,7:102])
colnames(datExpr)=colnames(dat1)[7:102]
rownames(datExpr)=dat1$Geneid
## Get the genes of interest
spikes <- rownames(datExpr)[grep("ERCC-", rownames(datExpr))]
genes <- rownames(datExpr)[!is.element(rownames(datExpr),spikes)]
## Also try with different numbers of factors
seqRUVg10 <- RUVg(seq, spikes, k=10)
## Write out the normalized expression data:
all.equal(as.character(rownames(normCounts(seqRUVg10))),as.character(dat1$Geneid))
#[1] TRUE
normCounts10=data.frame(dat1[,1:6],normCounts(seqRUVg10))
write.table(normCounts10,file="Normalized_read_counts_using_RUVg_ERCC_k10factors_outliers_removed.csv",col.names=TRUE,row.names=FALSE,sep=",")
# }}}
# run SampleNetwork
# {{{
setwd("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/RNAseq_non-normalized-hg19-ERCC92_trimmed_reads_stranded.csv")
datSample=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/SampleInformation_with_reads_hg19_ERCC92.csv")
datSample=data.frame(datSample, Group="All")
dim(dat1)
#[1] 26455   102
indexAll=c(7:102)
## Note: coloring by library concentration (22) to see whether there are any batch effects associated with concentration
dat1[,7:102]=apply(dat1[,7:102],2,as.numeric)
SampleNetwork(
datExprT=dat1,
method1="correlation",
impute1=FALSE,
subset1=NULL,
skip1=6,
indices1=list(indexAll),
subgroup1=22,
sampleinfo1=datSample,
samplelabels1=1,
grouplabels1=26,
fitmodels1=TRUE,
whichmodel1="univariate",
whichfit1="pc1",
btrait1=c(4,6,7,10,12,14,15,16,17,18,20,22,23,24),
trait1=NULL,
asfactors1=c(4,6,7,17,18,20,22,23,24),
projectname1="SF10711_mapped_to_hg19_and_ERCC92_stranded",
cexlabels=0.7,
normalize1=TRUE,
replacenegs1=FALSE,
verbose=TRUE
)
## After the first round remove <-2 which removed 6 samples. No more samples were removed after the first round. After quantile normalization there was a real outlyer which was the first sample. There was also a boarderline batch effect associated with RNA concentration. My bet is that this first sample was lower concentration so after removing it the batch effect will disappear. No batch correction was done after quantile normalization and a second round of SampleNetwork was run.
setwd("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/SF10711_mapped_to_hg19_and_ERCC92_stranded_SampleNetworks/All_10-53-46/SF10711_mapped_to_hg19_and_ERCC92_stranded_All_90_Qnorm.csv")
datSample=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/SampleInformation_with_reads_hg19_ERCC92.csv")
datSample=data.frame(datSample, Group="All")
dim(dat1)
#[1] 22969    96
datSample=datSample[is.element(datSample$SampleLabel,colnames(dat1)[7:96]),]
dim(datSample)
#[1] 90 26
indexAll=c(7:96)
## Note: coloring by library concentration (22) to see whether there are any batch effects associated with concentration
## Normalize set to false.
SampleNetwork(
datExprT=dat1,
method1="correlation",
impute1=FALSE,
subset1=NULL,
skip1=6,
indices1=list(indexAll),
subgroup1=22,
sampleinfo1=datSample,
samplelabels1=1,
grouplabels1=26,
fitmodels1=TRUE,
whichmodel1="univariate",
whichfit1="pc1",
btrait1=c(4,6,7,10,12,14,15,16,17,18,20,22,23,24),
trait1=NULL,
asfactors1=c(4,6,7,17,18,20,22,23,24),
projectname1="SF10711_mapped_to_hg19_and_ERCC92_stranded_round2",
cexlabels=0.7,
normalize1=FALSE,
replacenegs1=FALSE,
verbose=TRUE
)
## After the first round I used <-3 to remove sample 9-1-1. After removing this sample there were no batch effects. Potentially a minimal effect associated with sectioning plane but this is more likely to be biological related. Normalization or batch effects corrected for.
## }}} 
# run FindModules
# {{{
source("/fast-data/Shared/Code/KK/FindModules_0.90_KK.R")
setwd("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/FindModules")
dat1=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/RNAseq_data/Raw_data/fastq_files_101016/Analysis_in_R/RUVg_normalization/RUVg_normalized_k10Factors_SF10711_mapped_to_hg19_and_ERCC92_stranded_round3_SampleNetworks/All_01-01-26/RUVg_normalized_k10Factors_SF10711_mapped_to_hg19_and_ERCC92_stranded_round3_All_89_ComBat.csv")
FindModules(
projectname="SF10711_First_lane_RUVg_K=10",
expr=dat1,
geneinfo=c(1:6),
sampleindex=c(7:95),
samplegroups=NULL,
subset=NULL,
simMat=NULL,
saveSimMat=FALSE,
simType="Bicor",
beta=1,
overlapType="None",
TOtype="signed",
TOdenom="min",
MIestimator="mi.mm",
MIdisc="equalfreq",
signumType="rel",
iterate=TRUE,
signumvec=c(.9999,.999,.99,.98,.97,.96,.95,.94,.93,.92,.91,.90),
minsizevec=c(5,7,10,12,15),
signum=NULL,
minSize=NULL,
minMEcor=0.85,
ZNCcut=2,
calcSW=FALSE,
loadTree=FALSE,
writeKME=TRUE
)
# }}}
