# bcl2fastq
docker run --rm -v /tmp/data/input:/tmp/data/input:ro -v /tmp/data/output:/tmp/data/output umccr/bcl2fastq -R /tmp/data/input -o /tmp/data/output --sample-sheet /tmp/data/input/SampleSheet.csv --no-lane-splitting &> /tmp/data/output/bcl2fastq.log
sudo docker run --rm -v /home/patrick/190809_SF013_oldham_july_1/00.190802_NB502126_0210_AHVWG2BGXB:/home/patrick/190809_SF013_oldham_july_1/00.190802_NB502126_0210_AHVWG2BGXB:ro -v /home/patrick/output2:/home/patrick/output2 umccr/bcl2fastq -R /home/patrick/190809_SF013_oldham_july_1/00.190802_NB502126_0210_AHVWG2BGXB -o /home/patrick/output2 --minimum-trimmed-read-length 5 --mask-short-adapter-reads 5 --with-failed-reads --sample-sheet /home/patrick/190809_SF013_oldham_july_1/00.190802_NB502126_0210_AHVWG2BGXB/SampleSheet.csv --barcode-mismatches 1 --no-lane-splitting &> /home/patrick/output2/bcl2fastq.log
# correct umis with hamming distance of 2
for ea in *R1_001.fastq.gz; do 
	echo $ea
	sem -j 5 umi_tools whitelist --bc-pattern=CCCCCCCCCCCCNNNNNNNN --method=umis --knee-method=density --expect-cells=100 --error-correct-threshold=2 --ed-above-threshold=correct -L $ea.whitelist.log --stdin=$ea > $ea.whitelist ";" echo "done"
done
sem --wait
# extract umis and barcode
for ea in *R1_001.fastq.gz; do
	ea2=$( echo "$ea"  | sed 's/R1/R2/')
	sem -j 6 umi_tools extract --bc-pattern=CCCCCCCCCCCCNNNNNNNN --extract-method=string -L $ea.extract.log --error-correct-cell --whitelist $ea.whitelist  --filter-cell-barcode --read2-in=$ea2 --stdin=$ea --read2-out=$ea2.extract> /dev/null ";" echo "done"
done
# read trimming and adapter removal
for ea in *extract; do
java -jar /big-data2/Patrick/root_new/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 15 -trimlog $ea.log.trim  $ea $ea.trim ILLUMINACLIP:/big-data2/Patrick/root_new/bin/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:8 TRAILING:8 SLIDINGWINDOW:4:15 MINLEN:30 2>> $ea.trim.log
done
for ea in *extract; do
	TrimmomaticSE -threads 6 -trimlog $ea.log.trim  $ea $ea.trim ILLUMINACLIP:/usr/share/trimmomatic/NexteraPE-PE.fa:2:30:10 LEADING:8 TRAILING:8 SLIDINGWINDOW:4:15 MINLEN:30 2>> $ea.trim.log
done
# read alignment
for ea in *extract.trim; do 
	/big-data2/Patrick/root_new/bin/STAR-2.7.3a/bin/Linux_x86_64_static/STAR  --runThreadN 15 --genomeDir /big-data2/Patrick/GRCh38.p13.2019-02-28.gencode  --readFilesCommand cat --outSAMtype BAM SortedByCoordinate --readFilesIn $ea --outFilterMismatchNmax 999 --alignIntronMax 1000000 --outFilterMultimapNmax 20 --outSAMmultNmax 1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFileNamePrefix $ea.align
done
# read annotation
for ea in *.fastq.gz.extract.trim.alignAligned.sortedByCoord.out.bam; do
	/big-data2/Patrick/root_new/bin/subread-2.0.0-Linux-x86_64/bin/featureCounts -a /big-data2/Patrick/GRCh38.p13.2019-02-28.gencode/gencode.v33.chr_patch_hapl_scaff.annotation.ERCC.gtf -o $ea.assigned -R BAM  $ea -T 15 -g gene_id -M --primary -O -t gene
done
# sort and index bam file
for ea in *.featureCounts.bam; do
	/big-data2/Patrick/root_new/bin/samtools-1.10/samtools sort -@ 15 $ea >> $ea.sort
	/big-data2/Patrick/root_new/bin/samtools-1.10/samtools index -@ 15 $ea.sort
done
# get expression matrix {n.b. not used too slow}
for ea in *.extract.trim.alignAligned.sortedByCoord.out.bam.featureCounts.bam.sort; do 
	sem -j 6 umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I  $ea -S $ea.counts.tsv.gz >> log.umi.txt ";" echo "done"
done
