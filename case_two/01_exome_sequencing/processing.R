## Costello Lab Pipeline including Pindel output
dat1=read.delim("/Users/sheltons/Documents/Loglio_project/SF10711_astrocytoma/Exome_sequence_analysis/Exome_analysis_080516/Patient300.snvs.indels.filtered.txt")
dim(dat1)
#[1] 808 104
colnames(dat1)=gsub("X.gene","gene",colnames(dat1))
## Fix some of the gene symbols:
dat1$gene=gsub("(NM_014915:exon20:c.1986-2->TT)","",dat1$gene)
dat1$gene=gsub("(NM_018842:exon8:c.487-2->T)","",dat1$gene)
dat1$gene=gsub("C22orf29,","",dat1$gene)
dat1$gene=gsub("(NM_000722:exon9:c.659-2->T)","",dat1$gene)
dat1$gene=gsub("(NM_006587:exon6:c.618-2->T)","",dat1$gene)
dat1$gene=gsub(",CYP2D6","",dat1$gene)
dat1$gene=gsub("(NM_004447:exon21:c.2226-2->T)","",dat1$gene)
dat1$gene=gsub(",FRMD7","",dat1$gene)
dat1$gene=gsub("(NM_018072:exon46:c.6347-2->TT)","",dat1$gene)
dat1$gene=gsub(",KRT17","",dat1$gene)
dat1$gene=gsub("(NM_000081:exon19:c.5461-2->TT)","",dat1$gene)
dat1$gene=gsub(",MLH1","",dat1$gene)
dat1$gene=gsub("(NM_004145:exon16:c.2373+2->A,NM_001130065:exon16:c.2373+2->A)","",dat1$gene)
dat1$gene=gsub("(NM_012345:exon3:c.413-2->T)","",dat1$gene)
dat1$gene=gsub(",PDE1C","",dat1$gene)
dat1$gene=gsub("(NM_015617:exon3:c.50-2->T)","",dat1$gene)
dat1$gene=gsub("(NM_003821:exon8:c.1029+2->A)","",dat1$gene)
dat1$gene=gsub("(NM_145061:exon9:c.1120-2->TTT)","",dat1$gene)
dat1$gene=gsub(",TMEM161B","",dat1$gene)
dat1$gene=gsub(",UFL1","",dat1$gene)
dat1$gene=gsub("(NM_003400:exon9:c.591-2->T)","",dat1$gene)
dat1$gene=gsub("(NM_003821:exon8:c.1029+2->A)","",dat1$gene,fixed=TRUE)
dat1$gene=gsub("()","",dat1$gene,fixed=TRUE)
dat1$gene=gsub(",SELL","",dat1$gene)
dat1$gene=gsub("(NM_004145:exon16:c.2373+2->A,NM_001130065:exon16:c.2373+2->A)","",dat1$gene,fixed=TRUE)
## How many total mutations are there across all sample?
length(unique(dat1$gene))
#[1] 273
## There are still some silent mutations present in this list. Remove these as they are irrelevant
dat1=dat1[!is.na(!(dat1$type=="Silent"))&!(dat1$type=="Silent")|is.na(dat1$type),]
dim(dat1)
#[1] 627 104
length(unique(dat1$gene))
#[1] 225
## This reduced the number of mutations quite dramatically!
dat2=data.frame(table(as.character(dat1$gene)))
dat1Sub=dat1[is.element(dat1$gene,c("MUC17","TTN")),]
dat2a=table(as.character(dat1Sub$gene),dat1Sub$position)
dat2a
#100683876 100683944 179396304 179464170
#MUC17         1         1         0         0
#TTN           0         0         3         1
## The way to deal with these samples is to then sort by position and have multiple entries
datMT=data.frame(Var1=c("MUC17","MUC17","TTN","TTN"),Freq=c(1,1,3,1))
dat2=rbind(dat2[1:(grep("MUC17",dat2$Var1)-1),],datMT[1:2,],dat2[(grep("MUC17",dat2$Var1)+1):(grep("TTN",dat2$Var1)-1),],datMT[3:4,],dat2[(grep("TTN",dat2$Var1)+1):length(dat2[,1]),])
## At this point switch to the table that has just a single entry for each mutation:
dat1a=read.delim("/Users/sheltons/Documents/Loglio_project/SF10711_astrocytoma/Exome_sequence_analysis/Exome_analysis_080516/Patient300.snvs.indels.filtered.overlaps.txt")
dim(dat1a)
#[1] 277  87
colnames(dat1a)=gsub("X.gene","gene",colnames(dat1a))
dat1a$gene=gsub("(NM_014915:exon20:c.1986-2->TT)","",dat1a$gene)
dat1a$gene=gsub("(NM_018842:exon8:c.487-2->T)","",dat1a$gene)
dat1a$gene=gsub("C22orf29,","",dat1a$gene)
dat1a$gene=gsub("(NM_000722:exon9:c.659-2->T)","",dat1a$gene)
dat1a$gene=gsub("(NM_006587:exon6:c.618-2->T)","",dat1a$gene)
dat1a$gene=gsub(",CYP2D6","",dat1a$gene)
dat1a$gene=gsub("(NM_004447:exon21:c.2226-2->T)","",dat1a$gene)
dat1a$gene=gsub(",FRMD7","",dat1a$gene)
dat1a$gene=gsub("(NM_018072:exon46:c.6347-2->TT)","",dat1a$gene)
dat1a$gene=gsub(",KRT17","",dat1a$gene)
dat1a$gene=gsub("(NM_000081:exon19:c.5461-2->TT)","",dat1a$gene)
dat1a$gene=gsub(",MLH1","",dat1a$gene)
dat1a$gene=gsub("(NM_004145:exon16:c.2373+2->A,NM_001130065:exon16:c.2373+2->A)","",dat1a$gene)
dat1a$gene=gsub("(NM_012345:exon3:c.413-2->T)","",dat1a$gene)
dat1a$gene=gsub(",PDE1C","",dat1a$gene)
dat1a$gene=gsub("(NM_015617:exon3:c.50-2->T)","",dat1a$gene)
dat1a$gene=gsub("(NM_003821:exon8:c.1029+2->A)","",dat1a$gene)
dat1a$gene=gsub("(NM_145061:exon9:c.1120-2->TTT)","",dat1a$gene)
dat1a$gene=gsub(",TMEM161B","",dat1a$gene)
dat1a$gene=gsub(",UFL1","",dat1a$gene)
dat1a$gene=gsub("(NM_003400:exon9:c.591-2->T)","",dat1a$gene)
dat1a$gene=gsub("(NM_003821:exon8:c.1029+2->A)","",dat1a$gene,fixed=TRUE)
dat1a$gene=gsub("()","",dat1a$gene,fixed=TRUE)
dat1a$gene=gsub(",SELL","",dat1a$gene)
dat1a$gene=gsub("(NM_004145:exon16:c.2373+2->A,NM_001130065:exon16:c.2373+2->A)","",dat1a$gene,fixed=TRUE)
dat1a=dat1a[is.element(dat1a$gene,dat2$Var1),]
dim(dat1a)
#[1] 229  87
dat1a[duplicated(dat1a$gene),]
## Check these mutations.
dat1a[is.element(dat1a$gene,c("MUC17","HEATR1","TTN","EMR1")),]
## For these mutations there are multiple mutations for the same gene. For two of these (EMR1 and HEATR1) the second mutation is silent but for TTN and MUC17 there are two non-synonamous mutations so dat2 will need to be adjusted accordingly.
## Remove these from the data frame
dat1a=dat1a[!(dat1a$gene=="HEATR1"&dat1a$position==236749563),]
dim(dat1a)
#[1] 228  87
dat1a=dat1a[!(dat1a$gene=="EMR1"&dat1a$position==6924859),]
dim(dat1a)
#[1] 227  87
## Now sort these by gene and then by genomic position
dat1a=dat1a[order(as.character(dat1a$gene),as.numeric(dat1a$position)),]
dat2=dat2[order(as.character(dat2$Var1)),]
## Check that they are all in the correct order:
all.equal(as.character(dat1a$gene),as.character(dat2$Var1))
#[1] TRUE
## Combine these tables
dat1a=data.frame(dat1a[,1:10],no.samples.called.in=dat2$Freq,dat1a[,11:87])
## Now I would like to calculate the variant frequency for each of the mutations and add these to the table. It looks like there are a total of 13 exome samples that I will need to calculate the variant frequency for:
seq1=seq(19,74,4)
samples=gsub("_ref_reads","",colnames(dat1a)[seq1])
## Now create a data frame from these
datFreq=matrix(nrow=227,ncol=25)
colnames(datFreq)=c("gene","chr","position","ref_allele","alt_allele","nucleotide","protein","context","type","algorithm","no.samples.called.in",(paste(samples[1:14],".var.freq",sep="")))
## Now calculate variant frequencies for each sample in
for(i in 1:length(dat1a$gene)){
    mut1=dat1a[i,]
    if(mut1$algorithm=="Pindel"){
        datFreq[i,1]=as.character(mut1[1,1])
        datFreq[i,2]=as.character(mut1[1,2])
        datFreq[i,3]=as.numeric(mut1[1,3])
        datFreq[i,4]=as.character(mut1[1,4])
        datFreq[i,5]=as.character(mut1[1,5])
        datFreq[i,6]=as.character(mut1[1,6])
        datFreq[i,7]=as.character(mut1[1,7])
        datFreq[i,8]=as.character(mut1[1,8])
        datFreq[i,9]=as.character(mut1[1,9])
        datFreq[i,10]=as.character(mut1[1,10])
        datFreq[i,11]=as.numeric(mut1[1,11])
    }
    if(mut1$algorithm=="MuTect"){
        for(j in 1:length(samples)){
            if(j==2){
                datFreq[i,1]=as.character(mut1[1,1])
                datFreq[i,2]=as.character(mut1[1,2])
                datFreq[i,3]=as.numeric(mut1[1,3])
                datFreq[i,4]=as.character(mut1[1,4])
                datFreq[i,5]=as.character(mut1[1,5])
                datFreq[i,6]=as.character(mut1[1,6])
                datFreq[i,7]=as.character(mut1[1,7])
                datFreq[i,8]=as.character(mut1[1,8])
                datFreq[i,9]=as.character(mut1[1,9])
                datFreq[i,10]=as.character(mut1[1,10])
                datFreq[i,11]=as.numeric(mut1[1,11])
            }
            if(mut1[1,paste(samples[j],"_alt_Q20reads",sep="")]==0){
                datFreq[i,paste(samples[j],".var.freq",sep="")]=NA
            }else{
            datFreq[i,paste(samples[j],".var.freq",sep="")]=mut1[1,paste(samples[j],"_alt_Q20reads",sep="")]/(mut1[1,paste(samples[j],"_alt_Q20reads",sep="")]+mut1[1,paste(samples[j],"_ref_Q20reads",sep="")])
        }
    }
    }
}
## Now convert this into a data frame
datFreq=as.data.frame(datFreq)
## Write this file out to disk
write.table(datFreq,file="Summary_of_nonsynonamous_exonic_mutations_SF10711.csv",col.names=TRUE,row.names=FALSE,sep=",")
## Now I would like to create a bed file that I can use to subet down the bam files for each of these exome samples. I will then create an mpileup file for these and then use Varscan to try to get variant frequencies for each of these mutation (including the indels). To generate this bed file I will take 100bp upstream and downstream from the mutation site.
datFreq$position=as.character(datFreq$position)
datFreq$position=as.numeric(datFreq$position)
bedfile=data.frame(chrom=datFreq$chr,chromStart=(datFreq$position-100),chromEnd=(datFreq$position+100),name=paste(as.character(datFreq$gene),as.character(datFreq$type),sep="_"))
bedfile$name=gsub(" ","_",bedfile$name)
## Check that there are no other spaces:
bedfile$chrom=gsub(" ","_",bedfile$chrom)
bedfile$chromStart=gsub(" ","_",bedfile$chromStart)
bedfile$chromEnd=gsub(" ","_",bedfile$chromEnd)
## write this file out to disk. I will manually add a header line into it
write.table(bedfile,file="SF10711_exomeMutation.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
## How about running this to directly output a bam file.
## The -b flag means that the output is bam file. -L means it only outputs regions defined in the bed file.
samtools view -b -L SF10711_exomeMutation.bed.txt "$file" > ./Subsampled_bam_files/"${file/.bwa.realigned.rmDups.recal./_subsampled.}"
done
## Now sort the bam files (they should already by sorted but it is worth checking).
cd /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/Subsampled_bam_files
for file in *.bam
do
samtools sort "$file" ../Sorted_subsampled_bam_files/"${file/.bam/}".sorted
done
## Now index the sorted files:
cd /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/Sorted_subsampled_bam_files
for file in *.bam
do
samtools index "$file"
done
bedtools getfasta -fi hg19.fa -bed /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/Subsampled_reference_genome/SF10711_exomeMutation.bed.txt -fo subset_Hg19.fa
samtools faidx subset_Hg19.fa
##############
## Now make a mpileup file from the subset bam files. I will use a Q20 score since the Costello lab use this for their mutation calling.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/Sorted_subsampled_bam_files
mkdir ../mpileup
filelist=`ls -LR | grep .bam$`
filelist2=`echo $filelist | awk '{print}' ORS=' '`
filelist3=`echo $filelist2 | sed 's/*$//'`
samtools mpileup -A -d 999999 -Q 20 -f /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/Subsampled_reference_genome/hg19.fa $filelist3 >../mpileup/subset.SF10711.mpileup
mkdir ../SampleOrder
ls -LR | grep .bam$| sed 's/_subsampled.*//'| awk '{print}' ORS=',' | sed 's/*$//' >../SampleOrder/fileOrder.csv
## Now run Varscan to call the variants.
cd /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/mpileup
/fast-data/bin/jdk1.8.0_31/bin/java -jar /fast-data/bin/VarScan.v2.3.7.jar mpileup2cns /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/mpileup/subset.SF10711.mpileup --min-coverage 20 --min-var-freq 0.02 --p-value 0.20 --somatic-p-value 0.1 --validation 1 --variants >../VarScanOutput/Varscan_mpileup2cns_SF10711.txt
/fast-data/bin/jdk1.8.0_31/bin/java -jar /fast-data/bin/VarScan.v2.3.7.jar mpileup2cns /big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/mpileup/subset.SF10711.mpileup --min-coverage 20 --min-var-freq 0.08 --p-value 0.05 --somatic-p-value 0.1 --validation 1 --variants >../VarScanOutput/Varscan_mpileup2cns_SF10711-2.txt
## Now try to cross reference these lists with the Costello lab output.
setwd("/big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/VarScanOutput")
dat1=read.delim("/big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/VarScanOutput/Varscan_mpileup2cns_SF10711-2.txt")
refList=read.csv("/big-data2/Sam/Loglio_project/SF10711_astrocytoma/Exome_sequence_analysis/Exome_analysis_080516/Summary_of_nonsynonamous_exonic_mutations_SF10711.csv")
dim(dat1)
#[1] 329  11
dim(refList)
#[1] 227  25
## Use position to match the mutations
sum(is.element(dat1$Position,refList$position))
## It looks like 129 of the mutations have been matched by Varscan.
dat1Sub=dat1[is.element(dat1$Position,refList$position),]
refListSub=refList[is.element(refList$position,dat1Sub$Position),]
dim(dat1Sub)
#[1] 129  11
dim(refListSub)
#[1] 129  25
## Get them in the same order:
dat1Sub=dat1Sub[order(dat1Sub$Chrom,dat1Sub$Position),]
refListSub=refListSub[order(refListSub$chr,refListSub$position),]
all.equal(dat1Sub$Position,refListSub$position)
#[1] TRUE
all.equal(dat1Sub$Chrom,refListSub$chr)
#[1] TRUE
## Combine the tables together:
datCombi1=data.frame(gene=refListSub$gene,chr=refListSub$chr,postion=refListSub$position,Exome_ref_allele=refListSub$ref_allele,Exome_alt_allele=refListSub$alt_allele,Exome_nucleotide=refListSub$nucleotide,Exome_protein=refListSub$protein,context=refListSub$context,type=refListSub$type,Exome_algorthm=refListSub$algorithm,Exome_no.samples.called.in=refListSub$no.samples.called.in,Varscan.ref_allele=dat1Sub$Ref, Varscan.alt_allele=dat1Sub$Var, Varscan.no.samples.ref=dat1Sub$SamplesRef,Varscan.no.samples.het=dat1Sub$SamplesHet,Varscan.no.samples.hom=dat1Sub$SamplesHom,Varscan.no.samples.NC=dat1Sub$SamplesNC)
## To make the output easier to deal with in excell I will change + to ins and - to del.
datCombi1$Varscan.alt_allele=gsub("+","del",datCombi1$Varscan.alt_allele, fixed=TRUE)
datCombi1$Varscan.alt_allele=gsub("-","ins",datCombi1$Varscan.alt_allele, fixed=TRUE)
## Write out the data frame
write.table(datCombi1,file="Combined_Exome_Varscan_calls.csv",col.names=TRUE,row.names=FALSE,sep=",")
## Now determine the variant frequencies for each of the mutations
SampOrder=read.csv("/big-data2/Sam/Loglio_raw_data/SF10711/Exome_data/Recalibrated_bam_files_from_Tali/SampleOrder/fileOrder.csv",header=FALSE,nrow=1)
dim(SampOrder)
#[1]  1 14
SampOrder=t(SampOrder)
## This has split the data into 9 lists (1 per mutation) with 70 samples per list. I will deal with these 9 lists separately.
name<-paste("Mut",c(1:129),"df",sep="")
## I couldn't see ATRX, CTNND2, RECQL and WASL mutations.
prot=sapply(strsplit(as.character(refListSub$protein),split=",",fixed=TRUE),"[",1)
Muts<-c(paste(refListSub$gene,prot,sep="."))
SamplSplit=strsplit(as.character(dat1Sub[,11]),split=" ")
sampNames<-split(SampOrder, rep(1:ncol(SampOrder), each = nrow(SampOrder)))
for (i in 1:length(Muts)){
    Mut1Split=strsplit(SamplSplit[[i]],split=":")
    cons1=sapply(Mut1Split, "[", 1)
    cov1=sapply(Mut1Split, "[", 2)
    reads1a=sapply(Mut1Split, "[", 3)
    reads2a=sapply(Mut1Split, "[", 4)
    freq1=sapply(Mut1Split, "[", 5)
    freq1=gsub("%","",freq1)
    p.value1=sapply(Mut1Split, "[", 6)
    assign(name[i],data.frame(section.no=sampNames,cons=cons1,cov=cov1,reads1=reads1a,reads2=reads2a,p.cent.var.freq=freq1,p.val=p.value1))
}
## Rename the colnames of each data frame with the specific mutations from Muts:
mymuts <- lapply( name, get )
names(mymuts) <- Muts
for( j in length(mymuts) ) {
    colnames(mymuts[[j]])[c(2,3,4,5,6,7)] <- paste(Muts[j],colnames(mymuts[[j]])[c(2,3,4,5,6,7)],sep=".")
}
test<-do.call(cbind,mymuts)
SampLabel=as.character(test[,1])
test=test[,grep(".var.freq",colnames(test))]
test<-data.frame(lapply(test, as.character), stringsAsFactors=FALSE)
test<-data.frame(lapply(test, as.numeric), stringsAsFactors=FALSE)
test<-test/100
test=t(test)
colnames(test)=paste(SampLabel,".var.freq",sep="")
rownames(test)=gsub(".p.cent.var.freq","",as.character(rownames(test)),fixed=TRUE)
## Now combine these variant frequencies with the other data frame:
all.equal(as.character(datCombi1$gene),as.character(sapply(strsplit(rownames(test),split=".",fixed=TRUE),"[",1)))
[1] "1 string mismatch"
setdiff(as.character(datCombi1$gene),as.character(sapply(strsplit(rownames(test),split=".",fixed=TRUE),"[",1)))
#[1] "KRTAP5-1"
## For some reason this gene has been changed in test:
rownames(test)=gsub("KRTAP5.1","KRTAP5-1",rownames(test),fixed=TRUE)
all.equal(as.character(datCombi1$gene),as.character(sapply(strsplit(rownames(test),split=".",fixed=TRUE),"[",1)))
#[1] TRUE
datCombi1=data.frame(datCombi1,test)
write.table(datCombi1,file="Combined_Exome_Varscan_calls_with_variant_frequencies.csv",col.names=TRUE,row.names=FALSE,sep=",")
