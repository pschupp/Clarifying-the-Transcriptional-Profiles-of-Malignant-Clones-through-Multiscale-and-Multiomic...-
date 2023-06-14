## The raw exome data was processed through the our exome pipeline and the output from this pipeline was a list of somatic mutations that were identified by comparison to the exome data from the blood.
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data")
mut1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data/mutations.snvs.indels.csv")
dim(mut1)
# [1] 647  58
## I would like to get all of the mutations that are a real mutation by selecting the mutations that have "ourJudgment" column =yes.
# subselect real mutations and section name housekeeping
# {{{
mut1a = mut1[mut1[, 38] == "yes", ]
dim(mut1a)
# [1] 108  58
## There are 108 mutations in total that according to the Costello lab algorithm are genuine mutations.
length(unique(mut1a$gene))
# [1] 61
unique(mut1a$gene)
# [1] ABCC2    ACCS     ANKRD18B BRD7     C2       CALHM1   CSAG1    CTNNA3   FAM193B  GDPD2    GPR173
# [12] HIST1H4K IAPP     IDH1     IMMT     KRTAP5-3 KRTAP5-4 MUC16    MUC21    PASD1    PEX1     PHF8
# [23] PLA2G4D  SCN10A   SFTPA1   SIRT5    SUSD4    TBCC     TEC      TP53     TRIL     CRIPAK   CTSE
# [34] FAM60A   GPR126   NOL7     NTRK1    PABPC4L  PKD2     RECQL    RTN4     SCAPER   SLC20A2  SPEF2
# [45] WASL     XDH      APOBR    ATRX     C14orf39 CCDC168  CTNND2   CYSLTR2  DBC1     IL7R     MORC1
# [56] OR52N2   OR7C2    PCLO     SIGLEC5  SLAMF1   ZFYVE16
# 392 Levels: AAA1 ABCC2 ABI2 ABP1 ACACB ACCN5 ACCS ACSS2 ACTR3B ADCY8 AGAP11 AIM1L ... ZNF761
## There are 61 unique genes that are mutated in all three of the samples. Many of these mutations are likely to be the same across all of the samples. Some of these mutations are actually silent somatic mutations so they can be removed from the list:
dim(mut1a)
# [1] 108  58
mut1a = mut1a[!(mut1a$type == "Silent"), ]
dim(mut1a)
# [1] 89 58
## Also remove genes whose type is unknown. ANKRD18B is the only mutation with unknown type and I checked this muation and it is actually silent.
mut1a = mut1a[!(mut1a$type == "unknown"), ]
dim(mut1a)
# [1] 86 58
unique(mut1a$gene)
# [1] ABCC2    ACCS     BRD7     C2       CSAG1    CTNNA3   FAM193B  GDPD2
# [9] GPR173   HIST1H4K IAPP     IDH1     KRTAP5-3 MUC16    MUC21    PASD1
# [17] PEX1     PHF8     PLA2G4D  SUSD4    TBCC     TEC      TP53     TRIL
# [25] CRIPAK   FAM60A   GPR126   NOL7     NTRK1    PKD2     RECQL    RTN4
# [33] SCAPER   SLC20A2  SPEF2    WASL     ATRX     C14orf39 CCDC168  CTNND2
# [41] DBC1     IL7R     MORC1    OR52N2   OR7C2    PCLO     SLAMF1   ZFYVE16
## This leaves 48 non-synonamous somatic mutations.
length(unique(mut1a$gene))
# [1] 48
## This file still has the original sample numbers. I would like to amend this file to have the sequential renumbering. Read in the legend file that has both the old and new section numbers:
renum = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/relabelled section 1to69.csv")
dim(renum)
# [1] 69  2
colnames(renum)
# [1] "Original" "New"
## The current section numbers are 14, 39 and 69. The 69 section is actually mislabelled and is section #70. I will look these up in the section number reference table:
r14 = renum[renum$Original == 14, ]
r39 = renum[renum$Original == 39, ]
r70 = renum[renum$Original == 70, ]
renum = rbind(r14, r39, r70)
renum
# Original New
# 10       14  10
# 34       39  34
# 59       70  59
# Relabel the mutation list with the correct sequential number:
colnames(mut1a) = gsub("Glioma_14", "SF9495_10", colnames(mut1a))
colnames(mut1a) = gsub("Glioma_39", "SF9495_34", colnames(mut1a))
colnames(mut1a) = gsub("Glioma_69", "SF9495_59", colnames(mut1a))
## Also change the sample type:
mut1a$sample_type = gsub("Glioma_14", "SF9495_10", mut1a$sample_type)
mut1a$sample_type = gsub("Glioma_39", "SF9495_34", mut1a$sample_type)
mut1a$sample_type = gsub("Glioma_69", "SF9495_59", mut1a$sample_type)
## Also change the label "Glioma_Blood" to "SF9495_Blood":
mut1a$patient_ID = gsub("Glioma_Blood", "SF9495_Blood", mut1a$patient_ID)
## Write this file out for future reference.
write.table(mut1a, file = "non-synonamous somatic mutations.csv", col.names = TRUE, row.names = FALSE, sep = ",")
## Now I would like to narrow down the list to just those mutations that have a probe detected above background in at least one of the 69 tumor sections on the microarray:
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2013-140 sample probe profile.csv")
## Read in the original SampleInformation file that contains the missing samples (2 references and section 48)
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/SF9495_SampleInfo_for_SampleNetwork.csv")
dim(dat1)
# [1] 47202   302
dim(datSample)
# [1] 72 21
## I need to recalculate the CountDetPval and AvgDetPval to take account of the 3 samples that were removed during SampleNetwork (2 reference samples and section #48:
SampExcl = datSample[c(grep("SF9495_48", datSample$SampleLabel), grep("SF9495_REF_1", datSample$SampleLabel), grep("SF9495_REF_2", datSample$SampleLabel)), "Core_ID"]
## Now find these columns from the expression data and remove them:
grep(SampExcl[1], colnames(dat1))
# [1] 127 128 129 130
grep(SampExcl[2], colnames(dat1))
# [1] 51 52 53 54
grep(SampExcl[3], colnames(dat1))
# [1] 191 192 193 194
## Now remove these samples from the expression data:
dim(dat1)
# [1] 47202   302
dat1 = dat1[, -(grep(SampExcl[1], colnames(dat1)))]
dat1 = dat1[, -(grep(SampExcl[2], colnames(dat1)))]
dat1 = dat1[, -(grep(SampExcl[3], colnames(dat1)))]
dim(dat1)
# [1] 47202   290
# }}}
# process exome data based on detection parameters
# {{{
## Take the average signal for each of these samples
seq1 = seq(3, 275, 4)
## Calculate the average detection P value by taking the next column over from avg signal (seq1+1):
AvgDetPval = apply(dat1[, seq1 + 1], 1, mean)
## What is the distribution of mean detection P-values?
hist(AvgDetPval)
## Now we ask what fraction of transcripts are detected at P < .05 in at least one sample:
CountDetPval = apply(dat1[, seq1 + 1], 1, function(x) {
    length(x[x < .05])
})
## What is the distribution of the number of samples for which a given transcript was detected above background levels?
hist(CountDetPval, breaks = 50)
## Recreate the expression data with the CountDetPval and AvgDetPval calculated with the outliers removed:
dat2 = data.frame(dat1[, 1:2], dat1[, 286], AvgDetPval, CountDetPval, dat1[, seq1])
colnames(dat2)[2] = "Gene"
colnames(dat2)[3] = "RefSeq_ID"
colnames(dat2) = gsub(".AVG_Signal", "", colnames(dat2))
## Now I have recalculated the CountDetPval and AvgDetPval for the 69 samples I will replace these in the SampleNetwork processed expression data:
dat3 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")
dat2 = dat2[order(dat2$PROBE_ID), ]
dat3 = dat3[order(dat3$PROBE_ID), ]
all.equal(dat2$PROBE_ID, dat3$PROBE_ID)
# [1] TRUE
## Now replace the columns in dat3 with those in dat2:
dat3[, c(4:5)] = dat2[, c(4:5)]
## Check that there is no CountDetPvals greater than 69 in dat3 now it has been relabelled:
range(dat3[, 5])
# [1]  0 69
## It has been correctly relabelled so I will rewrite this table and overwrite the original since it contains incorrect information:
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor")
write.table(dat3, file = "Renumbered_SF9495_ALL_69_Qnorm.csv", col.names = TRUE, row.names = FALSE, sep = ",")
## Clear memory:
rm(list = ls())
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data")
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_SampleInfo_for_SampleNetwork.csv")
datSample$SampleLabel = paste("SF9495_", datSample$section_number, sep = "")
mut1a = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data/non-synonamous somatic mutations.csv")
## Check that the samples are in the same order between the sample information and expression files:
all.equal(as.character(datSample$SampleLabel), as.character(colnames(dat1)[7:75]))
# [1] TRUE
## Now I can start to ask what the detection p Values are for each of the mutations.
dat1a = dat1[is.element(dat1$Gene, intersect(dat1$Gene, mut1a$gene)), ]
## Of the 48 mutations, 74 probes were identified on the array.
dim(dat1a)
# [1] 74 75
## Now get the probes that are detected above background in at least one of the 69 tumor sections:
dat1b = dat1a[dat1a$CountDetPval >= 1, ]
dim(dat1b)
# [1] 62 75
## There are 44 mutated genes that are detected above background in at least 1 section.
length(unique(dat1b$Gene))
# [1] 44
## I will sort these by number of samples they are detected in and then by average detection p value.
dat1b = dat1b[order(-dat1b$CountDetPval, dat1b$AvgDetPval), ]
## Write out the table.
write.table(dat1b, file = "Probes_for_non-synonamous_mutations_Det_Above_background(1Sample).csv", row.names = FALSE, col.names = TRUE, sep = ",")
## It would also be useful to have the number of samples the mutation was identified in.
Det1 = data.frame(summary(mut1a$gene)[summary(mut1a$gene) > 0])
Det1 = data.frame(Gene = rownames(Det1), No.samples = Det1[, 1])
## Reduce this list to just the mutations that are detected above background in at least one of the tumor sections.
Det1a = Det1[is.element(Det1$Gene, intersect(Det1$Gene, dat1b$Gene)), ]
Det1a = Det1a[order(Det1a$Gene), ]
## Now I need to add in the number of samples the mutation was detected in:
dat1b = dat1b[order(dat1b$Gene), ]
Det2 = data.frame(summary(dat1b$Gene)[summary(dat1b$Gene) > 0])
Det2 = data.frame(Gene = rownames(Det2), No.probes = Det2[, 1])
Det2 = Det2[order(Det2[, 1]), ]
## Check that the mutation list (Det1a) are in the same order as the detection p-value data (Det2).
all.equal(as.character(Det1a[, 1]), as.character(Det2[, 1]))
# [1] TRUE
## Now replicate the numbers of samples by the numbers of probes detected above background.
int = c(rep.int(Det1a[1, 2], Det2[1, 2]))
for (n in 2:length(Det1a[, 1])) {
    int = c(int, rep.int(Det1a[n, 2], Det2[n, 2]))
}
int
# [1] 3 1 2 2 1 3 1 1 1 3 3 1 1 3 3 2 2 3 1 3 3 3 3 3 3 3 1 1 1 3 2 1 1 1 1 3 1 1
# [39] 1 1 1 1 3 2 2 2 1 1 1 1 1 1 1 1 3 3 3 3 3 1 1 1
length(int)
# [1] 62
## Now I will add this to the probe information
dat1b = data.frame(dat1b[, 1:5], Mut.no.samples = int)
## I will now reorder this first by the number of samples the mutation was identified in, then by the number of samples the gene was expressed above background in and then by Avg detection P value.
dat1b = dat1b[order(-dat1b$Mut.no.samples, -dat1b$CountDetPval, dat1b$AvgDetPval), ]
## Write this file over the original one:
## Now get the variant frequency for each of the mutations in each of the three exomes:
dat2c = dat1b
## Add blank columns to the data frame that I can put the allele variant frequency for each mutation into.
dat2c = data.frame(dat2c, Mut.Frac.10 = NA, Mut.Frac.34 = NA, Mut.Frac.59 = NA)
for (i in 1:length(dat2c$Gene)) {
    mut1aSub = mut1a[is.element(mut1a$gene, intersect(mut1a$gene, dat2c[i, 2])), ]
    if (sum(mut1aSub$sample_type == "SF9495_10") > 0) {
        dat2c[i, 7] = mut1aSub[mut1aSub$sample_type == "SF9495_10", 35]
    } else {
        dat2c[i, 7] = NA
    }
    if (sum(mut1aSub$sample_type == "SF9495_34") > 0) {
        dat2c[i, 8] = mut1aSub[mut1aSub$sample_type == "SF9495_34", 35]
    } else {
        dat2c[i, 8] = NA
    }
    if (sum(mut1aSub$sample_type == "SF9495_59") > 0) {
        dat2c[i, 9] = mut1aSub[mut1aSub$sample_type == "SF9495_59", 35]
    } else {
        dat2c[i, 9] = NA
    }
}
## Now order the data frame by the number of samples that the mutation was detected in, then by the number of probes that the genes expression was detected above background in then by AvgDetPval.
dat2c = dat2c[order(-dat2c$Mut.no.samples, -dat2c$CountDetPval, dat2c$AvgDetPval), ]
write.table(dat2c, file = "Probe_Summary_non-synonamous_mutations_Det_Above_background(in1Sample).csv", row.names = FALSE, col.names = TRUE, sep = ",")
## When there are duplicate probes it makes no sense to keep all of the probes so I will just take the probe that was detected above background in the most samples:
rm(list = ls())
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data/Probe_Summary_non-synonamous_mutations_Det_Above_background(in1Sample).csv")
dim(dat1)
# [1] 62  9
length(unique(dat1$Gene))
## [1] 44
# }}}
# focus on most likely mutations and incorporate sanger validation
# {{{
# Get a list of all the unique 44 genes that are detected above background in at least one section:
list = unique(dat1$Gene)
list = list[order(list)]
dat1 = dat1[order(dat1$Gene), ]
dat1a = data.frame(dat1[1, ])
dat1a[1, ] = NA
for (i in 1:length(list)) {
    get1 = grep(list[i], dat1$Gene)
    ## If there is more than 1 probe present on the array sort by CountDetPval and lowest AvgDetPval and take the highest Count and lowest AvgPval probe
    if (length(get1) < 1) {
        get2 = dat1[c(get1), ]
        ## Sort the probes by CountDetPval and then by AvgDetPval:
        get2 = get2[order(-get2$CountDetPval, get2$AvgDetPval), ]
        ## Take the probe with the highest CountDetPval and lowest AvgDetPval:
        get2 = get2[1, ]
        dat1a = rbind(dat1a, get2)
        ## If only one probe detected in the array for a gene just take that.
    } else {
        dat1a = rbind(dat1a, dat1[get1[1], ])
    }
}
## Remove the initial placeholder:
dat1a = dat1a[-1, ]
## This is the 44 mutations whose gene is detected above background in at least 1 section.
dim(dat1a)
# [1] 44  9
colnames(dat1a)
# [1] "PROBE_ID"       "Gene"           "RefSeq_ID"      "AvgDetPval"
# [5] "CountDetPval"   "Mut.no.samples" "Mut.Frac.10"    "Mut.Frac.34"
# [9] "Mut.Frac.59"
## Keep Gene, RefSeq_ID, CountDetPval, Mut.no.sample,Mut.Frac.10,Mut.Frac.34,Mut.Frac.59
dat1a = dat1a[, c(2, 3, 4, 5, 6, 7, 8, 9)]
## Rename the colnames
colnames(dat1a)[6] = "Var.freq.10"
colnames(dat1a)[7] = "Var.freq.34"
colnames(dat1a)[8] = "Var.freq.59"
## Order by CountDetPval then by Mut.no.samples
dat1a = dat1a[order(-dat1a$CountDetPval, -dat1a$Mut.no.samples), ]
## Write this out to disk
write.table(dat1a, file = "mutation_detected_above_background_summary.csv", col.names = TRUE, row.names = FALSE, sep = ",")
## Clear memory:
rm(list = ls())
## The next thing to do is to incorporate the amplicon Sanger sequencing data into the exome analysis. The majority of these 44 mutations were Sanger sequenced and a reasonably high proportion were found to also be present in the blood and were therefore false positives. I manually removed the mutations that were found to be false positives by Sanger sequencing and this file can be used to subset this data down. I also want to create a master table of mutation information based on this validated list of mutations:
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data")
## Read in the newly generated data frame with mutation and the number of samples it was called in as well as how many samples it was detected above background in.
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data/mutation_detected_above_background_summary.csv")
dim(dat1)
# [1] 44  8
datMut = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data/Complete list of unfiltered mutations.csv")
datMut = datMut[order(datMut$Gene), ]
dat1 = dat1[order(dat1$Gene), ]
all.equal(dat1$Gene, datMut$Gene)
# [1] TRUE
## Add the verified by Sanger column:
dat1 = data.frame(dat1, Verified.by.Sanger = datMut$Sanger.verified)
## Now load the output from the Costello pipeline restricted to somatic non-synonamous mutations:
dat2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Exome_data/non-synonamous somatic mutations.csv")
dim(dat2)
# [1] 86 58
## Restrict the list of confident mutations to the 33 which I am very confident about (called in exome data and look positive in IGV and/or Sanger validated).
dat2 = dat2[is.element(dat2$gene, dat1$Gene), ]
dim(dat2)
# [1] 82 58
## There are duplicate entries in this table because it adds a row for each of the 3 exomes it was called in.
length(unique(dat2$gene))
# [1] 44
## I now need to summarize the table by determining how many of the exomes each mutation was found in and which ones it was detected in.
dat2 = dat2[order(dat2$gene), ]
## Get the summary of the mutations. The data frame for the mutations is in alphabetical order. Since there are no cases where there are more than one mutation in the same gene I can use gene to summarize the data:
dat2a = as.data.frame(summary(dat2$gene)[summary(dat2$gene) > 0])
colnames(dat2a) = "No.Samples"
sum(dat2a$No.Samples == 3)
# [1] 17
sum(dat2a$No.Samples == 2)
# [1] 4
sum(dat2a$No.Samples == 1)
# [1] 23
17 + 4 + 23
# [1] 44
## Now I need to find out which of the exomes the mutation was called in. For the samples that are called in all 3 I don't need to worry about these. For the mutations called in 1 or 2 exomes I will determine which exomes the mutation was called.
dat2b = as.data.frame(summary(dat2$gene)[summary(dat2$gene) > 0 & summary(dat2$gene) < 3])
dat2b = data.frame(Gene = rownames(dat2b), No.Samples = dat2b[, 1])
calledin = c(1:length(dat2b[, 1]))
for (i in 1:length(dat2b[, 1])) {
    get1 = dat2[is.element(dat2$gene, dat2b[i, 1]), ]
    get2 = get1$sample_type
    get3 = gsub("SF9495_", "", get2)
    get3 = as.numeric(get3)
    if (length(get3) > 1) {
        get4 = paste(get3, collapse = ",")
        calledin[i] = as.character(get4)
    } else {
        calledin[i] = as.character(get3)
    }
}
dat2b = data.frame(dat2b, Called.In = calledin)
## Now combine this with the samples that were detected in all of the 3 exomes:
dat2c = as.data.frame(summary(dat2$gene)[summary(dat2$gene) > 2])
dat2c = data.frame(Gene = rownames(dat2c), No.Samples = dat2c[, 1], Called.In = paste(c(10, 34, 59), collapse = ","))
## Now combine the two data frames together:
dat2d = rbind(dat2b, dat2c)
dat2d = dat2d[order(as.character(dat2d$Gene)), ]
dim(dat2d)
# [1] 44  3
## Rename the columns:
colnames(dat2d) = c("Gene", "Mutated.No.Samples", "Mutation.Called.In")
## I also need to add what the actual mutation is into this table both the nucleotide and the protein change. There is an issue with the way that the data is displayed. When there are multiple transcripts for a gene there are multiple nucleotide and protein changes. This means that there are a number of different positions for the change depending on the transcript but the exome pipeline outputs all of them. For C2 and TP53 in one of the samples the order in which these variants are listed is different but they are the same variants just the orders are different. R sees these as different mutations but they are not.
## Change these manually:
# C2 first:
dat2[8:10, 1:6]
# gene contig position ref_allele alt_allele       nucleotide
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A
# 59   C2   chr6 31901429          G          A G485A,G116A,G89A
## Note the difference in order of the C2 nucleotide. Now I will change the order so they match the other 2:
dat2[10, "nucleotide"] = "G116A,G89A,G485A"
dat2[8:10, 1:6]
# gene contig position ref_allele alt_allele       nucleotide
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A
# 59   C2   chr6 31901429          G          A G116A,G89A,G485A
## Now for TP53:
dat2[75:77, 1:6]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 53 T38C,T434C,T38C,T434C,T38C,T434C,T434C
# 84 T38C,T434C,T434C,T38C,T38C,T434C,T434C
## Not that they are all in slightly different order so I will change them all to the first one:
dat2[75:77, "nucleotide"] = "T38C,T434C,T38C,T434C,T434C,T38C,T434C"
dat2[75:77, 1:6]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 53 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 84 T38C,T434C,T38C,T434C,T434C,T38C,T434C
## The same is true for the order of the protein change so I will correct that as well.
dat2[8:10, 1:7]
# gene contig position ref_allele alt_allele       nucleotide         protein
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 59   C2   chr6 31901429          G          A G116A,G89A,G485A R162Q,R39Q,R30Q
## Note the difference in order of the C2 protein. Now I will change the order so they match the other 2:
dat2[10, "protein"] = "R39Q,R30Q,R162Q"
dat2[8:10, 1:7]
# gene contig position ref_allele alt_allele       nucleotide         protein
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 59   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
## And for TP53
dat2[75:77, 1:7]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 53 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 84 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# protein
# 23 L13P,L145P,L13P,L145P,L145P,L13P,L145P
# 53 L13P,L145P,L13P,L145P,L13P,L145P,L145P
# 84 L13P,L145P,L145P,L13P,L13P,L145P,L145P
## Note that they are all in slightly different order so I will change them all to the first one:
dat2[75:77, "protein"] = "L13P,L145P,L13P,L145P,L145P,L13P,L145P"
dat2[75:77, 1:7]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 53 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 84 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# protein
# 23 L13P,L145P,L13P,L145P,L145P,L13P,L145P
# 53 L13P,L145P,L13P,L145P,L145P,L13P,L145P
# 84 L13P,L145P,L13P,L145P,L145P,L13P,L145P
## Since the order of the variants is also dependent on the accession numbers I need to change these as well.
dat2[8:10, c(1:7, 12)]
# gene contig position ref_allele alt_allele       nucleotide         protein
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 59   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# accession
# 4  NM_001178063,NM_001145903,NM_000063
# 27 NM_001178063,NM_001145903,NM_000063
# 59 NM_000063,NM_001178063,NM_001145903
dat2[9:10, "accession"] = "NM_001178063,NM_001145903,NM_000063"
dat2[8:10, c(1:7, 12)]
# gene contig position ref_allele alt_allele       nucleotide         protein
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 59   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# accession
# 4  NM_001178063,NM_001145903,NM_000063
# 27 NM_001178063,NM_001145903,NM_000063
# 59 NM_001178063,NM_001145903,NM_000063
## Now for TP53:
dat2[75:77, c(1:6, 12)]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 53 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 84 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# accession
# 23 NM_001126115,NM_001126114,NM_001126117,NM_001126113,NM_000546,NM_001126116,NM_001126112
# 53 NM_001126115,NM_001126114,NM_001126117,NM_001126113,NM_001126116,NM_000546,NM_001126112
# 84 NM_001126117,NM_001126113,NM_000546,NM_001126116,NM_001126115,NM_001126114,NM_001126112
dat2[76:77, "accession"] = "NM_001126115,NM_001126114,NM_001126117,NM_001126113,NM_000546,NM_001126116,NM_001126112"
dat2[75:77, c(1:6, 12)]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 53 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# 84 T38C,T434C,T38C,T434C,T434C,T38C,T434C
# accession
# 23 NM_001126115,NM_001126114,NM_001126117,NM_001126113,NM_000546,NM_001126116,NM_001126112
# 53 NM_001126115,NM_001126114,NM_001126117,NM_001126113,NM_000546,NM_001126116,NM_001126112
# 84 NM_001126115,NM_001126114,NM_001126117,NM_001126113,NM_000546,NM_001126116,NM_001126112
## The exonic position is also dependent on the accession so exon also needs to be altered.
dat2[8:10, c(1:7, 37)]
# gene contig position ref_allele alt_allele       nucleotide         protein
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 59   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# exon
# 4  2,2,4
# 27 2,2,4
# 59 4,2,2
dat2[9:10, "exon"] = "2,2,4"
dat2[8:10, c(1:7, 37)]
# gene contig position ref_allele alt_allele       nucleotide         protein
# 4    C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 27   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# 59   C2   chr6 31901429          G          A G116A,G89A,G485A R39Q,R30Q,R162Q
# exon
# 4  2,2,4
# 27 2,2,4
# 59 2,2,4
## Now for TP53:
dat2[75:77, c(1:6, 37)]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide          exon
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C 1,5,1,5,5,1,5
# 53 T38C,T434C,T38C,T434C,T434C,T38C,T434C 1,5,1,5,1,5,5
# 84 T38C,T434C,T38C,T434C,T434C,T38C,T434C 1,5,5,1,1,5,5
dat2[76:77, "exon"] = "1,5,1,5,5,1,5"
dat2[75:77, c(1:6, 37)]
# gene contig position ref_allele alt_allele
# 23 TP53  chr17  7578496          A          G
# 53 TP53  chr17  7578496          A          G
# 84 TP53  chr17  7578496          A          G
# nucleotide          exon
# 23 T38C,T434C,T38C,T434C,T434C,T38C,T434C 1,5,1,5,5,1,5
# 53 T38C,T434C,T38C,T434C,T434C,T38C,T434C 1,5,1,5,5,1,5
# 84 T38C,T434C,T38C,T434C,T434C,T38C,T434C 1,5,1,5,5,1,5
## I need to collapse dat2 to a single entry for each mutation. There are a couple of columns in the table that are unique to the sample in which it was called for example the variant frequency. Other than that the rest of the information in the table is duplicated. Since I am already including the variant frequency in the table anyway I can simply collapse the table but will use just the first few columns to collapse on:-
dat2 = dat2[!duplicated(dat2[, 1:7]), ]
dim(dat2)
# [1] 44 58
## Order the rows by gene alphabetical order:
dat2 = dat2[order(as.character(dat2$gene)), ]
all.equal(as.character(dat2$gene), as.character(dat2d$Gene))
# [1] TRUE
## Also get dat1 into the same order as dat2 so I can combine the two:
dat1 = dat1[order(as.character(dat1$Gene)), ]
all.equal(as.character(dat2$gene), as.character(dat1$Gene))
# [1] TRUE
## Looking at the known variant status, it doesn't make sense because it is saying that IDH1 R132H is novel when it is one of the most common glioma variants. For this reason I should leave it out of this table.
## Now add the mutation information to dat2d:
dat2d = data.frame(dat2d, Expression.AvgDetPval = dat1$AvgDetPval, Expression.CountDetPval = dat1$CountDetPval, contig = dat2$contig, position = dat2$position, accession = dat2$accession, mutation.context = dat2$context, exon = dat2$exon, mutation.type = dat2$type, algorithm = dat2$algorithm, nucleotide = dat2$nucleotide, protein = dat2$protein, Sanger.Verified = dat1$Verified.by.Sanger, COSMIC_mutation_frequency = dat2$COSMIC_mutation_frequency, COSMIC_mutation_within_3bp_frequency = dat2$COSMIC_mutation_within_3bp_frequency, COSMIC_gene_frequency = dat2$COSMIC_gene_frequency, KINASE = dat2$KINASE, dat1[, 6:8], dat2[, 42:58])
all.equal(as.character(dat2d$Gene), as.character(datMut$Gene))
dat2d = data.frame(dat2d, Amplicon.Sequenced = datMut$Amplicon.verified)
## Write out the comprehensive table:
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Tables_for_paper")
write.table(dat2d, file = "Complete_mutation_information.csv", col.names = TRUE, row.names = FALSE, sep = ",")
## It would also be nice to have a more cut down version of this table. Keep Gene, Mutated.No.Samples, Mutation.Called.in, Expression.CountDetPVal, contig, postion, mutation type, nucleotide, protein and Mut.Frac and SangerVerified:
dat2e = dat2d[, c(1, 6, 7, 13, 14, 8, 11, 3, 20:22, 5, 15, 40)]
## Rename the columns:
colnames(dat2e) = c("Gene", "Chromosome", "Position", "Nucleotide", "Protein", "Accession", "Mutation.type", "Mutation.distribution", "Var.freq.10", "Var.freq.34", "Var.freq.59", "No.sections.DABG", "Sanger.verified", "Amplicon.verified")
write.table(dat2e, file = "Cut_down_mutation_information.csv", col.names = TRUE, row.names = FALSE, sep = ",")
# }}} 
