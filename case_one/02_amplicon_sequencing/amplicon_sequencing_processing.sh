# first batch
# {{{
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Oldham_SS781
## Unzip the data from Oldham_SS781 and copy the fastq file to Second_run_fastq_files:
for f in *.gz; do
STEM=$(basename "${f}" .gz)
gunzip -c "${f}" > /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Second_run_fastq_files/"${STEM}"
done
## Next I need to create the reference in the appropriate format for the bwa algorithm. I manually generated the fastq file in textedit by copying and pasting the amplicon sequences together:
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/bwa-0.7.8
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Reference_genome_bwa
bwa index Glioma_Amplicon_Index.fasta
## Now align the reads to my custom reference sequence using bwa mem.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/bwa-0.7.8
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Second_run_fastq_files
mkdir ../Aligned_files_bwa/Sam_files
for R1 in *R1*.fastq
do
echo $R1
R2=`echo $R1 | sed 's/_R1_/_R2_/'`
newname=`echo $R1 | sed 's/_R1_001.fastq//'`
bwa mem /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Reference_genome_bwa/Glioma_Amplicon_Index.fasta $R1 $R2 > ../Aligned_files_bwa/Sam_files/$newname.sam
done
## Now convert the sam files to bam files
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sam_files
mkdir ../Bam_files
for file in *.sam
do
samtools view -b -S -o ../Bam_files/"${file/.sam/}".bam "$file"
done
## Next the bam files need to be sorted.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Bam_files
mkdir ../Sorted_bam_files
for file in *.bam
do
samtools sort "$file" ../Sorted_bam_files/"${file/.bam/}".sorted
done
## Now index the sorted bam files:
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
for file in *.bam
do
samtools index "$file"
done
## I will create a new folder to put the undetermined read in because I don't want these to be included in the variant frequency determination.
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
mkdir ../Undetermined_reads
for file in Undetermined*
do
mv $file ../Undetermined_reads
done
## Now use readcounts to count the reads support the reference or alternative allele.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
mkdir ../VarScanOutput/Varscan_readcounts_high_depth
filelist=`ls -LR | grep .bam$`
filelist2=`echo $filelist | awk '{print}' ORS=' '`
filelist3=`echo $filelist2 | sed 's/*$//'`
for file in $filelist3
do
samtools mpileup -Q 30 -d 999999999 -f /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Reference_genome_bwa/Glioma_Amplicon_Index.fasta "$file" |awk '($4>0)' | /fast-data/bin/jdk1.8.0_31/bin/java -jar /fast-data/bin/VarScan.v2.3.7.jar readcounts --variants-file /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/variants_list/variants.txt --min-base-qual 30 --min-coverage 0 --output-file ../VarScanOutput/Varscan_readcounts_high_depth/"${file/.bam/}"_readcounts.txt
done
# }}}
# second batch 
# {{{
## 11/23/15 The second batch of PCR amplicons has just been sequenced. I will process these in the same way that I processed the first batch.
## Unzip the data:
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/SS1463_4
## Unzip the data from SS1463_4 and copy the fastq file to Second_amplicon_batch_fastq_files:
for f in *.gz; do
STEM=$(basename "${f}" .gz)
gunzip -c "${f}" > /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Second_amplicon_batch_fastq_files/"${STEM}"
done
## For the first batch of amplicons the Gladstone core trimmed the adaptor sequences off the reads. For this second batch they left those on so I need to remove them before aligning the reads to the reference genome.
## Add trim_galore to the PATH.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
## Add trim_galore and FastQC to the PATH
PATH=$PATH:/fast-data/bin/trim_galore
PATH=$PATH:/fast-data/bin/FastQC
PATH=$PATH:/fast-data/bin
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Second_amplicon_batch_fastq_files
## Run trim_galore to remove the adaptor sequences. It automatically detects the adaptor type (nextera) and trims them off. Using the --paired argument deals with these as pairs and makes sure that there are reads of sufficient lengths for each pair otherwise the pair is discarded.
for R1 in *R1*.fastq
do
echo $R1
R2=`echo $R1 | sed 's/_R1_/_R2_/'`
trim_galore --paired $R1 $R2 --o ../Adaptor_trimmed_fastq_files
done
## Rename the files so they are in a similar format to the previous batch of amplicons:
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Adaptor_trimmed_fastq_files
## Rename read 1
for R1 in *R1*.fq
do
echo $R1
newname1=`echo $R1 | sed 's/_R1_001_val_1.fq/_R1_001.fastq/'`
cp $R1 ../Renamed_adaptor_trimmed_fastq_files/$newname1
done
## Rename read 2
for R1 in *R1*.fq
do
echo $R1
R2=`echo $R1 | sed 's/_R1_/_R2_/'`
R2=`echo $R2 | sed 's/val_1/val_2/'`
newname2=`echo $R2 | sed 's/_R2_001_val_2.fq/_R2_001.fastq/'`
cp $R2 ../Renamed_adaptor_trimmed_fastq_files/$newname2
done
## Next I need to create the reference in the appropriate format for the bwa algorithm. I manually generated the fastq file in textedit by copying and pasting the amplicon sequences together.:
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/bwa-0.7.8
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Reference_genome_bwa
bwa index Amplicons_fasta_second_run_only_reformated.index.fasta
## Now align the reads to my custom reference sequence using bwa mem.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/bwa-0.7.8
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Renamed_adaptor_trimmed_fastq_files
for R1 in *R1*.fastq
do
echo $R1
R2=`echo $R1 | sed 's/_R1_/_R2_/'`
newname=`echo $R1 | sed 's/_R1_001.fastq//'`
bwa mem /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Reference_genome_bwa/Amplicons_fasta_second_run_only_reformated.index.fasta $R1 $R2 > ../Aligned_files_bwa/Sam_files/$newname.sam
done
## Try running this on the cluster:
## This was run for each sample. Then using samtools I created the bam file.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Aligned_files_bwa/Sam_files
for file in *.sam
do
samtools view -b -S -o ../Bam_files/"${file/.sam/}".bam "$file"
done
## Sort the bam files on the cluster:
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Aligned_files_bwa/Bam_files
for file in *.bam
do
samtools sort "$file" ../Sorted_bam_files/"${file/.bam/}".sorted
done
## Now the sorted bam files need to be indexed. The index files need to be in the same directory as their matching bam file.:
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
for file in *.bam
do
samtools index "$file"
done
## I will create a new folder to put the undetermined read in because I don't want these to be included in the variant frequency determination.
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
mkdir ../Undetermined_reads
for file in Undetermined*
do
mv $file ../Undetermined_reads
done
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
mkdir ../VarScanOutput/Varscan_readcounts_high_depth
filelist=`ls -LR | grep .bam$`
filelist2=`echo $filelist | awk '{print}' ORS=' '`
for file in $filelist2
do
samtools mpileup -Q 30 -d 999999999 -f /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Reference_genome_bwa/Amplicons_fasta_second_run_only_reformated.index.fasta "$file" |awk '($4>0)' | /fast-data/bin/jdk1.8.0_31/bin/java -jar /fast-data/bin/VarScan.v2.3.7.jar readcounts --variants-file /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/Second_batch_amplicons/Variants_list/variants.txt --min-base-qual 30 --output-file ../VarScanOutput/Varscan_readcounts_high_depth/"${file/.bam/}"_readcounts.txt
done
# }}}
# combined batch one and two
# {{{
## Back in terminal. I was having problems generating the pileup files which I think is due to corrupted bam files. To overcome this I will reprocess this data.
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Oldham_SS781
## Unzip the data from Oldham_SS781 and copy the fastq file to Second_run_fastq_files:
for f in *.gz; do
STEM=$(basename "${f}" .gz)
gunzip -c "${f}" > /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Second_run_fastq_files/"${STEM}"
done
## Next I need to create the reference in the appropriate format for the bwa algorithm. I manually generated the fastq file in textedit by copying and pasting the amplicon sequences together:
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/bwa-0.7.8
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Reference_genome_bwa
bwa index Glioma_Amplicon_Index.fasta
## Now align the reads to my custom reference sequence using bwa mem.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/bwa-0.7.8
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Second_run_fastq_files
mkdir ../Aligned_files_bwa/Sam_files
for R1 in *R1*.fastq
do
echo $R1
R2=`echo $R1 | sed 's/_R1_/_R2_/'`
newname=`echo $R1 | sed 's/_R1_001.fastq//'`
bwa mem /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Reference_genome_bwa/Glioma_Amplicon_Index.fasta $R1 $R2 > ../Aligned_files_bwa/Sam_files/$newname.sam
done
## Now convert the sam files to bam files
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sam_files
mkdir ../Bam_files
for file in *.sam
do
samtools view -b -S -o ../Bam_files/"${file/.sam/}".bam "$file"
done
## Next the bam files need to be sorted.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Bam_files
mkdir ../Sorted_bam_files
for file in *.bam
do
samtools sort "$file" ../Sorted_bam_files/"${file/.bam/}".sorted
done
## Now index the sorted bam files:
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
for file in *.bam
do
samtools index "$file"
done
## I will create a new folder to put the undetermined read in because I don't want these to be included in the variant frequency determination.
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
mkdir ../Undetermined_reads
for file in Undetermined*
do
mv $file ../Undetermined_reads
done
## Now use readcounts to count the reads support the reference or alternative allele.
PATH=/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/dell/srvadmin/bin:/home/sheltons/bin
PATH=$PATH:/fast-data/bin/samtools-0.1.19
cd /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Aligned_files_bwa/Sorted_bam_files
mkdir ../VarScanOutput/Varscan_readcounts_high_depth
filelist=`ls -LR | grep .bam$`
filelist2=`echo $filelist | awk '{print}' ORS=' '`
filelist3=`echo $filelist2 | sed 's/*$//'`
for file in $filelist3
do
samtools mpileup -Q 30 -d 999999999 -f /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/Reference_genome_bwa/Glioma_Amplicon_Index.fasta "$file" |awk '($4>0)' | /fast-data/bin/jdk1.8.0_31/bin/java -jar /fast-data/bin/VarScan.v2.3.7.jar readcounts --variants-file /fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Amplicon_sequencing_data/First_batch_amplicons/variants_list/variants.txt --min-base-qual 30 --min-coverage 0 --output-file ../VarScanOutput/Varscan_readcounts_high_depth/"${file/.bam/}"_readcounts.txt
done
# }}}
