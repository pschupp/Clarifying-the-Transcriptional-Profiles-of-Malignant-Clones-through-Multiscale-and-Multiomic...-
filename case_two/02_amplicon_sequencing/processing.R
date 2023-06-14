# fastqc 
# {{{
for file in *.fastq
do
    fastqc "$file" --outdir ../Fastqc_prior_to_trimming
done
# }}}
# trim the reads 
# {{{
for file in *.fastq
do
(trim_galore -q 20 --length 20 "$file" --o ../Adaptor_trimmed_fastq_files) &
if (( $(wc -w <<<$(jobs -p)) % 10 == 0 )); then wait; fi # Limit to 10 concurrent subshells.
done
# }}}
# create custom reference and align
# {{{
cd /big-data2/Sam/Loglio_project/SF10711_astrocytoma/AmpSeq_data/AmpSeq_analysis/Reference_genome/
bwa index Combined_amplicons.fa
## Now align the reads to my custom reference sequence using bwa mem.
for R1 in *R1*.fq
do
    (newname=`echo $R1 | sed 's/_L001_R1_001_trimmed.fq//'`
    bwa mem /big-data2/Sam/Loglio_project/SF10711_astrocytoma/AmpSeq_data/AmpSeq_analysis/Reference_genome/Combined_amplicons.fa $R1  > ../Aligned_files_bwa_parallelized/Sam_files/$newname.sam) &
if (( $(wc -w <<<$(jobs -p)) % 10 == 0 )); then wait; fi # Limit to 10 concurrent subshells.
done
# }}}
# index and sort
# {{{
## Now convert the sam files to bam files
for file in *.sam
do
    (samtools view -b -S -o ../Bam_files/"${file/.sam/}".bam "$file") &
    if (( $(wc -w <<<$(jobs -p)) % 10 == 0 )); then wait; fi # Limit to 10 concurrent subshells.
done
wait
## Next the bam files need to be sorted.
for file in *.bam
do
    (samtools sort "$file" ../Sorted_bam_files/"${file/.bam/}".sorted) &
    if (( $(wc -w <<<$(jobs -p)) % 10 == 0 )); then wait; fi # Limit to 10 concurrent subshells.
done
wait
## Now index the sorted bam files:
for file in *.bam
do
    (samtools index "$file") &
    if (( $(wc -w <<<$(jobs -p)) % 10 == 0 )); then wait; fi # Limit to 10 concurrent subshells.
done
# }}}
# Get mutations
# {{{
#mpileup
for file in $filelist3
do
(samtools mpileup -Q 30 -d 999999999 -f /big-data2/Sam/Loglio_project/SF10711_astrocytoma/AmpSeq_data/AmpSeq_analysis/Reference_genome/Combined_amplicons.fa "$file" |awk '($4>0)' | /fast-data/bin/jdk1.8.0_31/bin/java -jar /fast-data/bin/VarScan.v2.3.7.jar readcounts --variants-file /big-data2/Sam/Loglio_project/SF10711_astrocytoma/AmpSeq_data/AmpSeq_analysis/Reference_genome/SF10711_variant_list.txt --min-base-qual 30 --min-coverage 0 --output-file ../VarScanOutput/Varscan_readcounts_high_depth/"${file/.bam/}"_readcounts.txt
) &
if (( $(wc -w <<<$(jobs -p)) % 10 == 0 )); then wait; fi # Limit to 10 concurrent subshells.
done
wait
# }}}
# process in R
# {{{
setwd("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/AmpSeq_data/AmpSeq_analysis/Aligned_files_bwa_parallelized/VarScanOutput/Varscan_readcounts_high_depth")
list1=list.files(getwd())
list1Name=paste(sapply(strsplit(list1,split="_",fixed=TRUE),"[",2),sapply(strsplit(list1,split="_",fixed=TRUE),"[",3),sapply(strsplit(list1,split="_",fixed=TRUE),"[",4),sep="_")
list1Name=data.frame(name=list1Name,Order=c(1:96))
NewOrder=sapply(strsplit(list1,split="_",fixed=TRUE),"[",2)
range(order(as.numeric(sapply(strsplit(list1[NewOrder=="9"],split="_",fixed=TRUE),"[",4)))+8)
#[1]  9 93
NewOrder[NewOrder=="10"]=94
NewOrder[NewOrder=="Primary"]=95
NewOrder[NewOrder=="Blood"]=96
NewOrderSub=as.numeric(sapply(strsplit(list1[NewOrder=="9"],split="_",fixed=TRUE),"[",4))
NewOrderSub=data.frame(section.no=NewOrderSub,order=c(1:85))
NewOrderSub=NewOrderSub[order(NewOrderSub$section.no),]
NewOrderSub=data.frame(NewOrderSub,New.order=c(9:93))
NewOrderSub=NewOrderSub[order(NewOrderSub$order),]
NewOrder[NewOrder=="9"]=NewOrderSub$New.order
NewOrder=as.numeric(NewOrder)
list1Name=data.frame(list1Name,NewOrder=NewOrder)
list1Name=list1Name[order(list1Name$NewOrder),]
## Now reorder the list1 to match list1Name
list1Name=list1Name[order(list1Name$Order),]
list1=list1[order(list1Name$NewOrder)]
list1Name=list1Name[order(list1Name$NewOrder),]
## I need to have a letter at the start of the name
list1Name$name=paste("SF10711",list1Name$name,sep="_")
all.equal(as.character(paste("SF10711",paste(sapply(strsplit(list1,split="_",fixed=TRUE),"[",2),sapply(strsplit(list1,split="_",fixed=TRUE),"[",3),sapply(strsplit(list1,split="_",fixed=TRUE),"[",4),sep="_"),sep="_")),as.character(list1Name$name))
#[1] TRUE
## Read all of these tables into R
for (i in 1:length(list1)){
    filetoread<-as.character(list1[i])
    ## To deal with the fact that the rows of the dataframe have different numbers of columns I will get the number of columns in each dataframe and create blank column names for them
    colnumb<-max(count.fields(filetoread, sep=""))
    assign(as.character(list1Name[i,1]),read.table(filetoread,sep="",col.names=paste("v",1:colnumb,sep="."),fill=T))
    assign(as.character(list1Name[i,1]),get(as.character(list1Name[i,1]))[2:dim(get(as.character(list1Name[i,1])))[1],1:13])
}
## I have a problem in that there are some samples that have 75 and some 73 mutations
getDim=matrix(ncol=2,nrow=96)
colnames(getDim)=c("Sample","No.Muts")
for(p in 1:length(list1Name[,1])){
    getDim[p,"Sample"]=list1Name[p,1]
    getDim[p,"No.Muts"]=dim(get(as.character(list1Name[p,1])))[1]
}
getDim=as.data.frame(getDim)
getDim[!getDim[,2]==74,]
#Sample No.Muts
#3  SF10711_3_1_19      73
#37 SF10711_9_1_46      75
#60 SF10711_9_2_81      73
#71 SF10711_9_2_95      73
## For the sake of time the easiest thing to do is to change all of these to 74 mutations by adding empty rows or deleting the extra row.
## Which of these are different
setdiff(as.character(SF10711_9_1_4[,1]),as.character(SF10711_3_1_19[,1]))
#[1] "MIS18A_A629G::chr21:33641346-33641496"
setdiff(as.character(SF10711_9_1_46[,1]),as.character(SF10711_9_1_4[,1]))
#[1] "FAM90A1_G368A::chr12:8376034-8376184"
setdiff(as.character(SF10711_9_1_4[,1]),as.character(SF10711_9_2_81[,1]))
#[1] "MIS18A_A629G::chr21:33641346-33641496"
setdiff(as.character(SF10711_9_1_4[,1]),as.character(SF10711_9_2_95[,1]))
#[1] "MIS18A_A629G::chr21:33641346-33641496"
## For 3-1-19, 9-2-81 and 9-2-95 I will add MIS18A. For 9-1-46 I will remove FAM90A1
SF10711_9_1_46=SF10711_9_1_46[!SF10711_9_1_46[,1]=="FAM90A1_G368A::chr12:8376034-8376184",]
dim(SF10711_9_1_46)
#[1] 74 13
SF10711_3_1_19=rbind(SF10711_3_1_19,SF10711_9_1_4[SF10711_9_1_4[,1]=="MIS18A_A629G::chr21:33641346-33641496",])
SF10711_3_1_19=SF10711_3_1_19[order(as.character(SF10711_3_1_19[,1])),]
SF10711_9_2_81=rbind(SF10711_9_2_81,SF10711_9_1_4[SF10711_9_1_4[,1]=="MIS18A_A629G::chr21:33641346-33641496",])
SF10711_9_2_81=SF10711_9_2_81[order(as.character(SF10711_9_2_81[,1])),]
SF10711_9_2_95=rbind(SF10711_9_2_95,SF10711_9_1_4[SF10711_9_1_4[,1]=="MIS18A_A629G::chr21:33641346-33641496",])
SF10711_9_2_95=SF10711_9_2_95[order(as.character(SF10711_9_2_95[,1])),]
## Check that the mutations are all in the same order
for (m in 1:length(list1)){
    if(m==1){
        MutOrder=cbind(as.character(get(as.character(list1Name[m,1]))[,1]),as.character(get(as.character(list1Name[m+1,1]))[,1]))
    }
    if(m>=3){
        MutOrder=cbind(MutOrder,as.character(get(as.character(list1Name[m,1]))[,1]))}
}
dim(MutOrder)
#[1] 74 96
## Check that all of the mutation in the first column match the order of the mutations in the next 69 columns.
## The mutations are in the correct order in the first two samples.
all.equal(as.character(MutOrder[,1]),as.character(MutOrder[,2]))
#[1] TRUE
## Now check the order in the remaining samples by using apply to match all of the columns to the first and count identical columns.
sum(apply(MutOrder[,2:96],2,identical, MutOrder[,1])==TRUE)
#[1] 95
## All of the mutations are in the same order and the data frames are the same dimension. There are 25 mutations and 70 samples. I want to make a data frame with samples (sections as rows) and columns as mutations. For each mutation I need reads for ref, reads for alt, total qX depth and frequency
## Create a blank matrix to populate with the required information
## Create mutation labels
MutLabel=sapply(strsplit(MutOrder[,1],split=":",fixed=TRUE),"[",1)
nmutations=length(MutLabel)
mat1=matrix(nrow=96,ncol=((nmutations*8)+1))
colnumber1=seq(2,dim(mat1)[2],8)
colnumber2=colnumber1+1
colnumber3=colnumber2+1
colnumber4=colnumber3+1
colnumber5=colnumber4+1
colnumber6=colnumber5+1
colnumber7=colnumber6+1
colnumber8=colnumber7+1
for (j in 1:dim(list1Name)[1]){
    if(j==1){
        mat1[,1]=list1Name[,1]
        colnames1=paste(rep(MutLabel,times=1,each=8),c("Q30.depth","Q30.reads.ref","Q30.reads.alt","freq","SE","LCI","UCI","conf.int"),sep="_")
        colnames(mat1)=c("sample",colnames1)
    }
    for(k in 1:nmutations){
        getstats=get(as.character(list1Name[j,1]))[k,c(5,6,8)]
        getstats[,2]<-sapply(strsplit(as.character(getstats[,2]),split=":",fixed=TRUE),"[",2)
        getstats<-as.numeric(t(getstats))
        ## Calculate the frequency based on reads supporting reference vs reads supporting ref vs alt
        frequency<-as.numeric(getstats[3])/(as.numeric(getstats[2])+as.numeric(getstats[3]))
        ## Calculate the SE of a percentage
        SE<-sqrt(frequency*(1-frequency)/getstats[1])
        ## Calculate the LCI, UCI and conf.int
        LCI=frequency-1.96*SE
        UCI=frequency+1.96*SE
        conf.int=UCI-LCI
        getstats=c(getstats,frequency,SE,LCI,UCI,conf.int)
        mat1[j,c(colnumber1[k],colnumber2[k],colnumber3[k],colnumber4[k],colnumber5[k],colnumber6[k],colnumber7[k],colnumber8[k])]<-getstats
    }
}
mat1=as.data.frame(mat1)
## write out this out to disk
setwd("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/AmpSeq_data/AmpSeq_analysis/Aligned_files_bwa_parallelized/Mutation_variant_frequency_data")
write.table(mat1,file="SF10711_variant_frequencies.csv",col.names=TRUE,row.names=FALSE,sep=",")
# }}}
