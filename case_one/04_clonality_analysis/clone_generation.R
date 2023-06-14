# Pipeline steps
# PyClone - cITUP
# [x] Pileup
# [x] FACETS
# [x] PyClone - 6 clusters
# [x] cITUP
# [x] Trees - waiting on PyCloneVI generated trees!
# [x] Chosing best
#
# PyClone-VI - cITUP
# [x] Pileup
# [x] FACETS
# [x] PyClone-VI - 8 clusters
# [x] cITUP - citup crashes with too many features. Currently ~400 features over 8 clusters seems to be the limit. But this is ok because after pyclone we input cellular prevalence and not vaf and the cellular prevalences are all extremly simliar anyway.  - use output_midi
# [x] Trees
# [x] Chosing best
# Mimicking the workflow used on the case SF9495, we will also try to derive a clonal tree from SF9495 and the related clonal prevalances across the three sections. In this case, it will be more difficult to identify cocorrelates, but perhaps by identifying key mutations in exonic sequences, we can try to then trace the VAF of mutations in the RNA-seq data and then identify the relevant module.
# GENERATING MPILEUP FILES
# {{{
# cd /mnt/bdata/patrick/SF9495/exome.seq/
# nice -n -18 samtools mpileup -B -q 1 -f /home/shared/hg_align_db/GRCh37.38/GRCh37.primary_assembly.genome.fa bamfiles/Glioma_Blood.bwa.realigned.rmDups.recal.bam | varscan mpileup2snp --min-freq-for-hom 0.95 --min-var-freq 0.05 --p-value 0.1 --output-vcf 1  > vcf_file/Glioma_Blood.bwa.realigned.rmDups.recal.bam.vcf
# Only SNPs will be reported
# Min coverage:   8
# Min reads2:     2
# Min var freq:   0.05
# Min avg qual:   15
# P-value thresh: 0.1
# Reading input from STDIN
# for ea in bamfiles/Glioma_[0-9]*bam; do
# 	echo $ea
# 	name=$(sed 's/.*\///g' <<< $ea)
# 	sem -j 4 nice -n -18 snp-pileup --count-orphans --max-depth=10000 --gzip --progress --min-map-quality=1  --min-base-quality=20 --min-read-counts 8,0 vcf_file/Glioma_Blood.bwa.realigned.rmDups.recal.bam.vcf pileup/$name.pileup bamfiles/Glioma_Blood.bwa.realigned.rmDups.recal.bam $ea
# done
# Calculating SNP count...done.
# Finished in 827.716820 seconds.========================================] 100 %
# double free or corruption (out)
# Calculating SNP count...done.
# Finished in 850.909731 seconds.========================================] 100 %
# double free or corruption (out)
# Calculating SNP count...done.
# Finished in 878.848829 seconds.========================================] 100 %
# double free or corruption (out)
source("/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R")
facets_run("/mnt/bdata/patrick/SF9495/exome.seq/", 50)
# }}}
# evidence of trunkal and clonal mutations
# {{{
# TRUNKAL
# 2p/q del at termini
# 4p del TEC # remove
# 7 amp PCLO
# 17p del  TP53 relative to rest of cluster
# CONAL
# 2p/q del at termini in 14 and 39
# 10p amp in 14 and 39
# 12p amp in 69 IAPP
# 12q del in 14
# }}}
# FOR PYCLONE, USING AMPSEQ SNPS + METHYLATION CNVS CALIBRATED BY EXOME CNVS
# {{{
source("/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R")
library("stringr")
amp = read.table("/mnt/bdata/patrick/SF9495/ampseq/Combined_variant_frequencies_run1_and_run2_with_readCounts.csv", sep = ",", header = T)
amp = amp[-(nrow(amp)), ]
colnames(amp) = gsub("DBC1", "BRINP1", colnames(amp))
colnames(amp) = gsub("GPR126", "ADGRG6", colnames(amp))
colnames(amp) = gsub("HIST1H4K", "H4C12", colnames(amp))
pur14 = read.csv("/mnt/bdata/patrick/SF9495/exome_seq/Glioma_14.stats.csv")
pur39 = read.csv("/mnt/bdata/patrick/SF9495/exome_seq/Glioma_39.stats.csv")
pur69 = read.csv("/mnt/bdata/patrick/SF9495/exome_seq/Glioma_69.stats.csv")
# copy number derived from FACETS
cnv14 = read.csv("/mnt/bdata/patrick/SF9495/exome_seq/facets/Glioma_14.cval.250.csv")
cnv39 = read.csv("/mnt/bdata/patrick/SF9495/exome_seq/facets/Glioma_39.cval.250.csv")
cnv69 = read.csv("/mnt/bdata/patrick/SF9495/exome_seq/facets/Glioma_69.cval.250.csv")
cnv14$chrom = gsub(23, "X", cnv14$chrom)
cnv39$chrom = gsub(23, "X", cnv39$chrom)
cnv69$chrom = gsub(23, "X", cnv69$chrom)
# unique cnv parsing
temp1 = cnv14[which(cnv14$Minor.Copy.Number != 1 | cnv14$Major.Copy.Number != 1), ]
cnvs = data.frame(temp1, sec = rep("s14", nrow(temp1)))
temp1 = cnv39[which(cnv39$Minor.Copy.Number != 1 | cnv39$Major.Copy.Number != 1), ]
cnvs = rbind(cnvs, data.frame(temp1, sec = rep("s39", nrow(temp1))))
temp1 = cnv69[which(cnv69$Minor.Copy.Number != 1 | cnv69$Major.Copy.Number != 1), ]
cnvs = rbind(cnvs, data.frame(temp1, sec = rep("s69", nrow(temp1))))
cnvs = cnvs[order(cnvs$chrom, cnvs$start), ]
cnvs = rbind(cnvs, rep(0, ncol(cnvs)))
cnvs_u = data.frame(matrix(nrow = 1, ncol = ncol(cnvs)))
colnames(cnvs_u) = colnames(cnvs)
cnvs_u = cnvs_u[order(as.numeric(cnvs_u$chrom), cnvs_u$start), ]
# cnvs_u is a dataframe that lists all the unique cnvs which we believe to exist across the sections
j = 1
cnvCol = list()
for (i in seq(1, (nrow(cnvs) - 1))) {
    if (cnvs$chrom[i] == cnvs$chrom[i + 1]) {
        if (cnvs$start[i] >= (cnvs$start[i + 1] * .9) & cnvs$start[i] <= (cnvs$start[i + 1] * 1.1)) {
            if (cnvs$end[i] >= (cnvs$end[i + 1] * .9) & cnvs$end[i] <= (cnvs$end[i + 1] * 1.1)) {
                if (length(cnvCol) < j) {
                    cnvCol[[j]] = c("NA")
                }
                cnvCol[[j]] = c(cnvCol[[j]], cnvs$Cellular.Fraction[i])
                next
            }
        }
    }
    cnvs_u[j, ] = cnvs[i, ]
    if (length(cnvCol) < j) {
        cnvCol[[j]] = c("NA")
    }
    if (length(cnvCol) == 0) {
        next
    }
    cnvCol[[j]] = c(cnvCol[[j]], cnvs$Cellular.Fraction[i])
    names(cnvCol[[j]]) = c("NA", cnvs$sec[(i - (length(cnvCol[[j]]) - 2)):i])
    j = j + 1
}
cnvs_u = cnvs_u[order(as.numeric(cnvs_u$chrom), cnvs_u$start), ]
cnvCol = lapply(cnvCol, function(x) x[-1])
names(cnvCol) = paste(cnvs_u$chrom, cnvs_u$start, sep = "_")
cnvCol = lapply(cnvCol, function(x) x[order(names(x))])
cnvColDf = data.frame(matrix(nrow = 3, ncol = length(cnvCol)))
for (i in seq(1, length(cnvCol))) {
    a = cnvCol[[i]]
    if (length(grep("s69", names(a))) == 0) {
        a = c(a, 0)
        names(a)[length(a)] = "s69"
    }
    if (length(grep("s39", names(a))) == 0) {
        a = c(0, a)
        names(a)[c(1)] = c("s39")
    }
    if (length(grep("s14", names(a))) == 0) {
        a = c(0, a)
        names(a)[1] = "s14"
    }
    a = a[order(names(a))]
    cnvColDf[, i] = a
}
cnvColDf = apply(cnvColDf, 2, as.numeric)
rownames(cnvColDf) = c("14", "39", "69")
colnames(cnvColDf) = names(cnvCol)
# cnvColDf is a dataframe that lists the CNVs across sections derived by FACETS
# read in methylation data
cnvMeth = fread("~/bdata/SF9495/methylation_array/LogR_ratio_tumor_vs_normal_brain.csv")
cnvMethOut = data.frame(matrix(nrow = 1, ncol = ncol(cnvMeth)))
for (i in seq(1, nrow(cnvs_u))) {
    temp = cnvMeth[which(cnvMeth$chrom == cnvs_u$chrom[i] & cnvMeth$pos > cnvs_u$start[i] & cnvMeth$pos < cnvs_u$end[i]), ]
    temp2 = apply(temp[, -c(1, 2)], 2, median)
    cnvMethOut[i, ] = c(temp[1, c(1, 2)], temp2)
}
colnames(cnvMethOut) = colnames(cnvMeth)
rownames(cnvMethOut) = paste(cnvMethOut$chrom, cnvMethOut$pos, sep = "_")
cnvMethOut = cnvMethOut[, -c(1, 2)]
relabel = read.csv("~/bdata/SF9495/methylation_array/relabelled section 1to69.csv")
renameIndex = match(relabel$New, as.numeric(gsub(".*_", "", colnames(cnvMethOut))))
colnames(cnvMethOut) = relabel$Original[renameIndex[!is.na(renameIndex)]]
cnvMethOut = t(cnvMethOut)
# Ampseq is missing section 37
cnvMethOut = cnvMethOut[-which(rownames(cnvMethOut) == 37), ]
# convert methylation data cellular frequencies of CN
cnvMethOutCor = data.frame(matrix(nrow = nrow(cnvMethOut), ncol = ncol(cnvMethOut)))
for (i in seq(1, ncol(cnvMethOut))) {
    methWork = cnvMethOut[which(rownames(cnvMethOut) %in% rownames(cnvColDf)), i]
    cnvMethOutCor[, i] = cnvMethOut[, i] * (mean(cnvColDf[cnvColDf[, i] != 0, i] / methWork[cnvColDf[, i] != 0]))
}
cnvMethOutCor[cnvMethOutCor < 0] = 0
rownames(cnvMethOutCor) = rownames(cnvMethOut)
colnames(cnvMethOutCor) = c("chr2p_del", "chr2q_del", "chr4p_del", "chr7_amp", "chr10p_amp", "chr12p_amp", "chr12q_del", "chr17p_del")
# read in GTF data
gtf = fread("/home/shared/hg_align_db/GRCh37.38/gencode.v38lift37.annotation.gene_only.gtf", sep = "\t")
gtf = gtf[, -c(2, 3, 6, 8)]
colnames(gtf) = c("chr", "start", "end", "strand", "gene")
gtf$gene = gsub('";.*', "", gsub('.*gene_name "', "", gtf$gene))
genes = gsub("\\..*", "", colnames(amp)[grep("depth", colnames(amp))])
gtf = gtf[match(genes, gtf$gene), ]
gtf$chr = gsub("chr", "", gtf$chr)
gtf = gtf[order(gtf$chr, gtf$start), ]
gtf$mean = apply(gtf[, c(2, 3)], 1, mean)
# assign copy number to gtf output
gtf = cnAssign(cnv14, "s14")
gtf = cnAssign(cnv39, "s39")
gtf = cnAssign(cnv69, "s69")
# track distance of ampseq section to nearest exome section
sections = as.numeric(gsub(".*_", "", amp$sample))
sectionsDist = data.frame(abs(sections - 14), abs(sections - 39), abs(sections - 69))
sectionsMin = apply(sectionsDist, 1, which.min)
setwd("/mnt/bdata/patrick/SF9495/clonality/pyclone/inputs")
amp = amp[which(amp$sample %in% rownames(cnvMethOutCor)), ]
names = paste("sf", str_pad(gsub(".*_", "", amp$sample), width = 3, pad = 0, side = "left"), sep = "")
for (i in seq(1, nrow(amp))) {
    sample = amp$sample[i]
    mut = gsub("_Q.*", "", colnames(amp)[grep("depth", colnames(amp))])
    mut2 = gsub(".*_", "", mut)
    mut = gsub("\\..*", "", mut)
    mut = paste(mut, mut2, sep = "_")
    ref = amp[i, grep("\\.ref$", colnames(amp))]
    alt = amp[i, grep("\\.alt$", colnames(amp))]
    ind1 = as.numeric(6 + 3 * sectionsMin[i])
    ind2 = as.numeric(5 + 3 * sectionsMin[i])
    out = data.frame(mutation_id = mut, ref_counts = as.numeric(ref), var_counts = as.numeric(alt), normal_cn = rep(2, length(mut)), minor_cn = gtf[match(gsub("_.*", "", mut), gtf$gene), ..ind1], major_cn = gtf[match(gsub("_.*", "", mut), gtf$gene), ..ind2])
    colnames(out) = c("mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn")
    out$minor_cn[which(is.na(out$minor_cn))] = 0
    out$major_cn[which(is.na(out$major_cn))] = 2
    pseudosnps = data.frame(mutation_id = colnames(cnvMethOutCor), ref_counts = as.numeric(round(mean(out$ref_counts) * (1 - cnvMethOutCor[i, ]))), var_counts = as.numeric(round(mean(out$ref_counts) * (cnvMethOutCor[i, ]))), normal_cn = rep(2, ncol(cnvMethOutCor)), minor_cn = rep(0, ncol(cnvMethOutCor)), major_cn = rep(1, ncol(cnvMethOutCor)))
    out = rbind(out, pseudosnps)
    write.table(out, file = paste(names[i], "tsv", sep = "."), sep = "\t", row.names = F)
}
setwd("/mnt/bdata/patrick/SF9495/clonality/pyclone/inputs")
out_com = paste("/mnt/bdata/patrick/SF9495/clonality/pyclone/inputs/", names, ".tsv", sep = "", collapse = ",")
out_com = paste("PyClone run_analysis_pipeline --in_files {", out_com, "} ", "--tumour_contents {", paste(amp$TP53.137.TtoC_freq, collapse = ","), "} ", "--working_dir /mnt/bdata/patrick/SF9495/clonality/pyclone/outputs", sep = "")
setwd("/mnt/bdata/patrick/SF9495/clonality/pyclone/")
write.table(out_com, file = "pyclone_command.txt", quote = F, row.names = F, col.names = F)
# }}}
# MAKE SURE YOU ARE NOT USING X-FORWARDING
# {{{
# conda activate python2
# chmod +x pyclone_command.txt
# ./pyclone_command.txt
# conda deactivate
# }}}
# convert output to cITUP format
# {{{
cluster = fread("/mnt/bdata/patrick/SF9495/clonality/pyclone/outputs/tables/cluster.tsv")
loci = fread("/mnt/bdata/patrick/SF9495/clonality/pyclone/outputs/tables/loci.tsv")
loci = loci[order(loci$sample_id), ]
freq = reshape(loci[, c(1, 2, 4)], timevar = "sample_id", idvar = "mutation_id", direction = "wide")
clusters = aggregate(loci$cluster_id, by = list(loci$mutation_id), FUN = mean)
colnames(clusters) = c("mutation_id", "cluster")
clusters = clusters[order(clusters$cluster), ]
freq = freq[match(clusters$mutation_id, freq$mutation_id), ]
# remove singletons
singletons = clusters$mutation_id[which(clusters$cluster %in% names(table(clusters$cluster))[which(table(clusters$cluster) == 1)])]
singletons = singletons[c(3)] # chose which singletons to exclude
if (length(singletons) > 0) {
    freq = freq[!(freq$mutation_id == singletons), ]
    clusters = clusters[!(clusters$mutation_id == singletons), ]
}
write.table(clusters$cluster, quote = F, row.names = F, col.names = F, sep = "\t", file = "/mnt/bdata/patrick/SF9495/clonality/pyclone_cITUP/inputs/clusters.txt")
write.table(freq[, -1], quote = F, row.names = F, col.names = F, sep = "\t", file = "/mnt/bdata/patrick/SF9495/clonality/pyclone_cITUP/inputs/freq.txt")
# }}}
# run cITUP data
# {{{
# ALWAYS RUN WITH MIN/MAX NODES CLONES +1 BECAUSE BASE NODE IS PARENTAL!
# cd /mnt/bdata/patrick/SF9495/clonality/pyclone_cITUP/outputs
# conda activate pyclone
# run_citup_qip.py --maxjobs 5 --nocleanup --min_nodes 7 --max_nodes 7 ../inputs/freq.txt ../inputs/clusters.txt results.h5
# conda deactivate
# }}}
# Plot tree from cITUP output
# {{{
source("/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R")
citupToTree(citDir = "/mnt/bdata/patrick/SF9495/clonality/pyclone_cITUP/outputs", snpDir = "/mnt/bdata/patrick/SF9495/clonality/pyclone/outputs/tables", outDir = "/mnt/bdata/patrick/SF9495/clonality/pyclone_cITUP/trees")
clusterW = cast(cluster, sample_id ~ cluster_id, value = "mean")
clusterW[, 6] = clusterW[, 6] - clusterW[, 2]
clusterW[, 2] = clusterW[, 2] - apply(clusterW[, 3:5], 1, sum)
clusterW[, 3:5] = clusterW[, 3:5] * .75
colnames(clusterW) = c("sample_id", 2, 3, 5, 1, 4)
clusterM = melt(clusterW)
clusterM$cluster_id = as.factor(clusterM$cluster_id)
clusterM[clusterM < 0] = 0
pdf("try.pdf")
print(ggplot(clusterM, aes(x = sample_id, y = value, fill = cluster_id)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)))
dev.off()
# }}}
# FOR PYCLONE-VI, USING EXOME DATA
# {{{
# convert pilup locations to genes and protein changes
# --offline
# cd /mnt/bdata/patrick/SF9495/exome_seq/pileup/
# for(ea in list.files()[grep('gz', list.files())]){
# 	pileup=fread(ea)
# 	pileup=pileup[-which(pileup$Chromosome %in% c('chrM', 'chrY','chrX')),]
# 	pileupDepth=apply(pileup[,c('File2R','File2A')], 1, sum)
# 	pileup=pileup[(pileupDepth>200),]
# 	pileOut=data.frame(chr=gsub('chr', '', pileup$Chromosome), start=pileup$Position, stop=pileup$Position, allele=paste(pileup$Ref, pileup$Alt, sep='/'))
# 	fwrite(pileOut, row.names=F, col.names=F, sep='\t', quote=F, file=paste0(gsub('gz', '', ea), 'vep.gz'), compress='gzip')
# }
# for ea in *vep.gz; do
# 	eaN=${ea/.*/}
# 	/opt/vep_ensembl/ensembl-vep/vep	-i $ea -o $eaN.vep.txt --cache --assembly GRCh37 --offline --dir_cache /home/shared/vep_cache/ \
# 										--fasta /home/shared/vep_cache/fasta/GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa --hgvs --hgvsg --shift_hgvs 1 \
# 										--symbol --refseq --fork 1 --pick --tab \
# 										--fields  "Uploaded_variation,SYMBOL,STRAND,Feature,HGVSc,HGVSp,HGVSg"
# done
# }}}
# Reading in pileup files
# {{{
cnAssignNumr = function(mut, cnv) {
    major_cn = minor_cn = rep(NA, nrow(mut))
    for (i in seq(1, nrow(cnv14))) {
        ind = which(mut$Chromosome == cnv$chrom[i] & mut$Position >= cnv$start[i] & mut$Position <= cnv$end[i])
        if (length(ind) == 0) {
            next
        }
        major_cn[ind] = rep(cnv$Major.Copy.Number[i], length(ind))
        minor_cn[ind] = rep(cnv$Minor.Copy.Number[i], length(ind))
    }
    major_cn[is.na(major_cn)] = 1
    minor_cn[is.na(minor_cn)] = 1
    mut = data.frame(mut, major_cn = major_cn, minor_cn = minor_cn)
    return(mut)
}
mut14 = read.table("/mnt/bdata/patrick/SF9495/exome.seq/pileup/Glioma_14.bwa.realigned.rmDups.recal.bam.pileup.gz", sep = ",", header = T, colClasses = c("factor", "numeric", "character", "character", "NULL", "NULL", "NULL", "NULL", "numeric", "numeric", "NULL", "NULL"))
mut39 = read.table("/mnt/bdata/patrick/SF9495/exome.seq/pileup/Glioma_39.bwa.realigned.rmDups.recal.bam.pileup.gz", sep = ",", header = T, colClasses = c("factor", "numeric", "character", "character", "NULL", "NULL", "NULL", "NULL", "numeric", "numeric", "NULL", "NULL"))
mut69 = read.table("/mnt/bdata/patrick/SF9495/exome.seq/pileup/Glioma_69.bwa.realigned.rmDups.recal.bam.pileup.gz", sep = ",", header = T, colClasses = c("factor", "numeric", "character", "character", "NULL", "NULL", "NULL", "NULL", "numeric", "numeric", "NULL", "NULL"))
totalReads = data.frame(chr = gsub("chr", "", mut14$Chromosome), s14 = apply(mut14[, 5:6], 1, sum), s39 = apply(mut39[, 5:6], 1, sum), s69 = apply(mut69[, 5:6], 1, sum))
totalReads$mean = apply(totalReads[, 2:4], 1, mean)
mut14 = mut14[which(totalReads$mean > 200), ] # for fewer mutations
mut39 = mut39[which(totalReads$mean > 200), ] # for fewer mutations
mut69 = mut69[which(totalReads$mean > 200), ] # for fewer mutations
mut14$Chromosome = gsub("chr", "", mut14$Chromosome)
mut39$Chromosome = gsub("chr", "", mut39$Chromosome)
mut69$Chromosome = gsub("chr", "", mut69$Chromosome)
cnv14 = read.csv("/mnt/bdata/patrick/SF9495/exome.seq/facets/Glioma_14.cval.250.csv")
cnv39 = read.csv("/mnt/bdata/patrick/SF9495/exome.seq/facets/Glioma_39.cval.250.csv")
cnv69 = read.csv("/mnt/bdata/patrick/SF9495/exome.seq/facets/Glioma_69.cval.250.csv")
mut14 = cnAssignNumr(mut = mut14, cnv = cnv14)
mut39 = cnAssignNumr(mut = mut39, cnv = cnv39)
mut69 = cnAssignNumr(mut = mut69, cnv = cnv69)
out = rbind(mut14, mut39, mut69)
out$mutation_id = paste(out$Chromosome, out$Position, sep = ":")
out$sample_id = rep(c("mut14", "mut39", "mut69"), ea = nrow(mut14))
out$normal_cn = rep(2, nrow(out))
out = out[, c(10, 11, 5, 6, 12, 8, 9, 7)]
colnames(out) = c("mutation_id", "sample_id", "ref_counts", "alt_counts", "normal_cn", "major_cn", "minor_cn", "tumour_content")
out = out[-grep("^M|^Y|^X", out$mutation_id), ]
write.table(out, file = "/mnt/bdata/patrick/SF9495/clonality/pyclone-vi/sf9495_input_mini.tsv", quote = F, sep = "\t", row.names = F)
# }}}
# RUNNING PYCLONE-VI
# {{{
# cd /mnt/bdata/patrick/SF9495/clonality/pyclone-vi/outputs
# activate pyclone-vi
# pyclone-vi fit -i sf9495_input.tsv -o sf9495_output.h5 -c 40 -d beta-binomial -r 10
# pyclone-vi write-results-file -i sf9495_output.h5 -o sf9495_output.tsv
# conda deactivate
# }}}
# read in vep files  which list genes and consequences from snp pileup files
# {{{
s14 = fread("/mnt/bdata/patrick/SF9495/exome_seq/pileup/Glioma_14.vep.txt", skip = "#Uploaded_variation")
s39 = fread("/mnt/bdata/patrick/SF9495/exome_seq/pileup/Glioma_39.vep.txt", skip = "#Uploaded_variation")
s69 = fread("/mnt/bdata/patrick/SF9495/exome_seq/pileup/Glioma_69.vep.txt", skip = "#Uploaded_variation")
# condense all unique mutations
snpS = list()
snpS[[1]] = which(!(s39$"#Uploaded_variation" %in% s14$"#Uploaded_variation"))
snpS[[2]] = which(!(s69$"#Uploaded_variation" %in% s14$"#Uploaded_variation"))
ind = which(colnames(s14) %in% c("#Uploaded_variation", "SYMBOL", "HGVSp", "HGVSg"))
snpLeg = rbind(s14[, ..ind], s39[snpS[[1]], ..ind], s69[snpS[[2]], ..ind])
colnames(snpLeg) = c("mutation_id", "gene", "hgvsp", "hgvsg")
ind = which(snpLeg$hgvsp == "-")
snpLeg$names = snpLeg$hgvsp
snpLeg$names[ind] = snpLeg$hgvsg[ind]
snpLeg$names = gsub(".*\\.", "", snpLeg$names)
snpLeg$names = paste(snpLeg$gene, snpLeg$names, sep = ":")
snpLeg = snpLeg[, c("mutation_id", "names")]
snpLeg$mutation_id = gsub("_", ":", gsub("_[A-z].*", "", snpLeg$mutation_id))
# read in pyclone data
pyclone = fread("/mnt/bdata/patrick/SF9495/clonality/pyclone-vi/outputs/sf9495_output_mini.tsv")
# change pyclone names to useful gene names
ind = which(pyclone$mutation_id %in% snpLeg$mutation_id)
pyclone$mutation_id[ind] = snpLeg$names[match(pyclone$mutation_id[ind], snpLeg$mutation_id)]
pyclone$mutation_id = make.names(pyclone$mutation_id)
ids = paste0(pyclone$mutation_id, pyclone$sample_id)
pyclone$mutation_id[duplicated(ids)] = paste0(pyclone$mutation_id[duplicated(ids)], ".1")
toRemove = grep("^X\\.\\.", pyclone$mutation_id)
toRemove = c(toRemove, grep("^LOC", pyclone$mutation_id))
toRemove = c(toRemove, grep("^X[0-9]", pyclone$mutation_id))
pyclone = pyclone[-toRemove, ]
# }}}
# MAX 400 FEATURES across all clones
# {{{
pycloneMini = data.frame(matrix(nrow = 1, ncol = ncol(pyclone)))
colnames(pycloneMini) = colnames(pyclone)
for (ea in unique(pyclone$cluster_id)) {
    temp = pyclone[pyclone$cluster_id == ea, ]
    if (nrow(temp) > 150) {
        rand = sample(seq(1, nrow(temp), 3), 50)
        temp = temp[sort(c(rand, rand + 1, rand + 2)), ]
    }
    pycloneMini = rbind(pycloneMini, temp)
}
# for midi, cluster 6 has 28000 members, reduce this to the max of the second most populous cluster 4, which has 4746
freq = reshape(pycloneMini[, c(1, 2, 4)], idvar = "mutation_id", timevar = "sample_id", direction = "wide")
clusters = aggregate(pycloneMini$cluster_id, by = list(pycloneMini$mutation_id), FUN = mean)
colnames(clusters) = c("mutation_id", "cluster")
clusters = clusters[order(clusters$cluster), ]
freq = freq[match(clusters$mutation_id, freq$mutation_id), ]
# remove singletons
singletons = clusters$mutation_id[which(clusters$cluster %in% names(table(clusters$cluster))[which(table(clusters$cluster) == 1)])]
singletons = singletons[c()] # chose which singletons to exclude
if (length(singletons) > 0) {
    freq = freq[!(freq$mutation_id == singletons), ]
    clusters = clusters[!(clusters$mutation_id == singletons), ]
}
write.table(clusters$cluster, quote = F, row.names = F, col.names = F, sep = "\t", file = "/mnt/bdata/patrick/SF9495/clonality/pyclone-vi_cITUP/inputs/clusters_midi.txt")
write.table(freq[, -1], quote = F, row.names = F, col.names = F, sep = "\t", file = "/mnt/bdata/patrick/SF9495/clonality/pyclone-vi_cITUP/inputs/freq_midi.txt")
# }}}
# run cITUP data
# {{{
# ALWAYS RUN WITH MIN/MAX NODES CLONES +1 BECAUSE BASE NODE IS PARENTAL!
# cd /mnt/bdata/patrick/SF9495/clonality/pyclone-vi_cITUP/outputs
# conda activate pyclone
# run_citup_qip.py --maxjobs 15 --min_nodes 8 --max_nodes 8 --nocleanup ../inputs/freq_midi.txt ../inputs/clusters_midi.txt results.h5
# conda deactivate
# }}}
# Plot tree from cITUP output
# {{{
source("/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R")
citupToTree(citDir = "/mnt/bdata/patrick/SF9495/clonality/pyclone_cITUP/outputs", snpDir = "/mnt/bdata/patrick/SF9495/clonality/pyclone/outputs/tables", outDir = "/mnt/bdata/patrick/SF9495/clonality/pyclone_cITUP/trees")
# }}}
# Make clone barplot
# {{{
clusterW = cast(cluster, sample_id ~ cluster_id, value = "mean")
clusterW[, 6] = clusterW[, 6] - clusterW[, 2]
clusterW[, 2] = clusterW[, 2] - apply(clusterW[, 3:5], 1, sum)
clusterW[, 3:5] = clusterW[, 3:5] * .75
colnames(clusterW) = c("sample_id", 2, 3, 5, 1, 4)
clusterM = melt(clusterW)
clusterM$cluster_id = as.factor(clusterM$cluster_id)
clusterM[clusterM < 0] = 0
pdf("try.pdf")
print(ggplot(clusterM, aes(x = sample_id, y = value, fill = cluster_id)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)))
dev.off()
# }}}
