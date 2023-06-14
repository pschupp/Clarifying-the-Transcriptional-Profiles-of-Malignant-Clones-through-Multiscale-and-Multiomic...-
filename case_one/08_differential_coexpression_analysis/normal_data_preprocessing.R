## I prepared raw sample information during the preparation of the samples. This includes RIN scores from the bioanalyzer, conc ng/ul, OD260/280 etc. I took this information and prepared it for use in SampleNetwork.

################################ Preprocessing the normal 4 datasets from the ACC and EC from 2 individuals ####################

## The ACC1 first:

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1")

## Note: added array info and position info from core spreadsheet to sample information file.  Reordered by Array_position column (e.g. 6303256015_A).

datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_SampleInfo/ACC1_SampleInfo_for_SampleNetwork.csv")
datSample = datSample[, -18]

## To translate the continuous first randomization number (QC_number) to a QC_batch:

QC_batch = ceiling((as.integer(datSample$QC_number) / 12))
datSample = data.frame(datSample, QC_batch)
## Write the amended table to the working directory.
write.table(datSample, file = "ACC1_SampleInfo_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)

## Note: open "2011-199 sample probe profile.txt" in Excel and save it as a .csv file
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2011-199 sample probe profile.csv")
dim(dat1)
# [1] 47231   398

## We will select the expression values ("AVG_Signal") columns and summarize the Detection P-values:
## We will also retain probe metadata.

seq1 = seq(3, 383, 4)

AvgDetPval = apply(dat1[, seq1 + 1], 1, mean)

## What is the distribution of mean detection P-values?

hist(AvgDetPval)

## So ~20K/47K transcripts show a mean detection P-value < .05 (i.e. are expressed significantly above background levels).

## Now we ask what fraction of transcripts are detected at P < .05 in at least one sample:

CountDetPval = apply(dat1[, seq1 + 1], 1, function(x) {
    length(x[x < .05])
})

## What is the distribution of the number of samples for which a given transcript was detected above background levels?

hist(CountDetPval, breaks = 50)

dat2 = data.frame(dat1[, 1:2], dat1[, 394], AvgDetPval, CountDetPval, dat1[, seq1])
colnames(dat2)[2] = "Gene"
colnames(dat2)[3] = "RefSeq_ID"
colnames(dat2) = gsub(".AVG_Signal", "", colnames(dat2))

## Now reorder datSample by CoreLabel:

datSample = datSample[order(datSample$CoreLabel), ]

## Check to see that expression data columns and sample information rows are in the same order:

all.equal(as.character(gsub("X", "", colnames(dat2)[6:101])), as.character(datSample$CoreLabel))
# [1] TRUE

## So all samples are in the same order.  Lets re-label the column headers in the expression data:

colnames(dat2)[6:101] = as.character(datSample$SampleLabel)

## Give the reference samples a section number so that they can be ordered correctly:

datSample[datSample$SampleLabel == "REF_1", "Section_number"] = 187
datSample[datSample$SampleLabel == "REF_2", "Section_number"] = 188

## Now get the columns in the correct order:

datSample = datSample[order(datSample$Section_number), ]
datSample = data.frame(datSample, new.order = c(1:96))
datSample = datSample[order(datSample$CoreLabel), ]

all.equal(as.character(colnames(dat2)[6:101]), as.character(datSample$SampleLabel))
# [1] TRUE

dat2Sub = dat2[, 6:101]
dat2Sub = dat2Sub[, order(datSample$new.order)]
datSample = datSample[order(datSample$Section_number), ]

all.equal(as.character(colnames(dat2Sub)), as.character(datSample$SampleLabel))
# [1] TRUE

dat2 = data.frame(dat2[, 1:5], dat2Sub)

## Reorder the expression data by PROBE_ID:

dat2 = dat2[order(dat2$PROBE_ID), ]

## Add the probe quality.

datAnnot = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1/ACC1_ExpressionData_for_SampleNetwork.csv")

datAnnot = datAnnot[order(datAnnot$PROBE_ID), ]

all.equal(datAnnot$PROBE_ID, dat2$PROBE_ID)
# [1] TRUE

dat2 = data.frame(dat2[, 1:2], ProbeQuality = datAnnot$ProbeQuality, dat2[, 3:101])

write.table(dat2, file = "ACC1_ExpressionData_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)
write.table(datSample, file = "ACC1_SampleInfo_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)

rm(list = ls())

## Run SampleNetwork (SampleNetwork run locally because it is easier to run interactively than on the cluster):

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1 = read.csv("ACC1_ExpressionData_for_SampleNetwork.csv")
datSample = read.csv("ACC1_SampleInfo_for_SampleNetwork.csv")

indexAll = c(7:102)

## Note: coloring by ArrayID to start (subgroup1=13), which may or may not be a significant batch effect.

SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 13,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 17,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 13, 14, 18, 11, 8, 5, 6),
    trait1 = NULL,
    asfactors1 = c(13, 14, 18),
    projectname1 = "ACC1",
    cexlabels = 0.7,
    normalize1 = TRUE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)

## Note: After the first round I removed the two reference samples by <-3. Second round one sample was removed using <-3 and one sample was removed in the third round using <-3 again. After normalization there was a significant batch effect associated with Array_ID which was corrected for using ComBat with no covariates.


## Preprocessing EC1.

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC1")
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_SampleInfo/EC1_SampleInfo_for_SampleNetwork.csv")

## The initial expression data sheet that I got back from the UCLA core had a mistake as some samples were mislabelled. They sent me a replacement called RS 357corrected.

dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_SampleInfo/RS 357corrected.csv")
colnames(dat1)

# [1] "pre.sort.RS.Index"                "Running.Set.Well.ID"
# [3] "SCGC.Project.ID"                  "PI.or.Lab"
# [5] "Project.Type"                     "Organism"
# [7] "Number.of.Samples"                "Sample.Type"
# [9] "Chip.Type"                        "External.ID"
# [11] "general.array"                    "General.Array.Stripe.entry.field"
# [13] "Date"                             "Primary.contact.name"
# [15] "Primary.contact.Email"            "Primary.Contact.Address"
# [17] "Phone"                            "Institute"
# [19] "Project.cost"                     "Account.numberbilling"
# [21] "Fund.numberbilling"               "cost.centerbilling"
# [23] "Fund.manager.emailbilling"        "LOCbilling"
# [25] "SUBbilling"                       "Objectbilling"
# [27] "projectbilling"                   "Running.Set"

colnames(datSample)

# [1] "Section_number" "QC_number"      "Core_number"    "ND_conc_ng.ul"  "ND_260.280"
# [6] "ND_260.230"     "BA_conc_ng.ul"  "RIN"            "Ratio_28s.18s"  "FirstScreen"
# [11] "MeanConc"       "PMI"            "QC_batch"

## Create a unique core number for dat1 by pasting the general array with General Array Stripe entry field.
corNum = paste(dat1[, 11], dat1[, 12], sep = "_")

## Add this to dat1:
dat1 = data.frame(dat1, CoreLabel = corNum)
## Remove the NA's at the end of dat1 which are just blank entries:
dim(dat1)
# [1] 100  29
dat1 = dat1[1:96, ]
dim(dat1)
# [1] 96 29
dat1 = data.frame(dat1, Order = c(1:96))

## Now order datSample by Core_number:
datSample = datSample[order(datSample$Core_number), ]

## Now move one of the reference samples to the start of the table:
refs = datSample[95:96, ]
datSample = datSample[1:94, ]
datSample = rbind(refs[1, ], datSample, refs[2, ])

## Now create SampleLabel:
SampleLabel = paste("EC", datSample[, 1], sep = "_")

SampleLabel[1] = "REF_1"
SampleLabel[96] = "REF_2"

datSample = data.frame(datSample, SampleLabel = SampleLabel)

## Relabel the reference samples in dat1:
dat1[, 10] = as.character(dat1[, 10])
dat1[1, 10] = NA
dat1[96, 10] = NA

## Check that the samples are in the same order:
all.equal(as.character(dat1[, 10]), as.character(datSample[, 1]))
# [1] TRUE

## Create QC batch.
QC_batch = ceiling((as.integer(datSample$QC_number) / 12))

## Now combine the two data tables:

datSample2 = data.frame(datSample[, 1:12], ArrayID = (dat1[1:96, 11]), ArrayPosition = (dat1[1:96, 12]), CoreLabel = dat1$CoreLabel, SampleLabel = datSample$SampleLabel, Group = rep("ALL", length = 96), QC_Batch = (QC_batch))

## Now order these by CoreLabel since the


## Write the new table as sample info for sample network:

write.table(datSample2, file = "EC1_SampleInfo_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)

rm(list = ls())

## Now reload the new SampleInfo file and read in the probe data. I had some issues writing this table out on my local machine. For some reason some rows were omitted from the output. I think this must be an R version number problem because I don't get this issue if I do it on the cluster with version 3.1.1 (my local machine has 3.1.2).

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC1")
datSample = read.csv("EC1_SampleInfo_for_SampleNetwork.csv")
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2011-242 sample probe profile.csv")
## Load the data from the ACC1:
dat2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2011-199 sample probe profile.csv")

## Order the data frames by PROBE_ID
dat1 = dat1[order(dat1[, 1]), ]
dat2 = dat2[order(dat2[, 1]), ]

dim(dat1)
# [1] 47231   391

dim(dat2)
# [1] 47231   398

all.equal(as.character(dat1[, 1]), as.character(dat2[, 1]))
# [1] TRUE
## The gene symbols are also the same
all.equal(as.character(dat1[, 2]), as.character(dat2[, 2]))
# [1] TRUE
## How about the REFSEQ_ID.
all.equal(as.character(dat1[, 387]), as.character(dat2[, 394]))
## [1] "12417 string mismatches"
## This looks to be the problem. For some reason there is no REFSEQ_ID instead there is something that looks like a refseq_ID as SEARCH_KEY. This is clearly different from that of the ACC1. Since the probes are exactly the same I will just use the REFSEQ_ID from the ACC1 expression data.

## We will select the expression values ("AVG_Signal") columns and summarize the Detection P-values By selecting the columns with data in them
## We will also retain probe metadata.

seq1 = seq(3, 386, 4)

AvgDetPval = apply(dat1[, seq1 + 1], 1, mean)

## What is the distribution of mean detection P-values?

hist(AvgDetPval)

## So ~20K/47K transcripts show a mean detection P-value < .05 (i.e. are expressed significantly above background levels).

## Now we ask what fraction of transcripts are detected at P < .05 in at least one sample:

CountDetPval = apply(dat1[, seq1 + 1], 1, function(x) {
    length(x[x < .05])
})

## What is the distribution of the number of samples for which a given transcript was detected above background levels?
hist(CountDetPval, breaks = 50)

## Since the probes are the same between the ACC1 data and EC1 (same platform) and the EC1 data lacks a refseq ID I will use the one from the ACC1 (dat2 column 394).
dat3 = data.frame(dat1[, 1:2], dat2[, 394], AvgDetPval, CountDetPval, dat1[, seq1])
colnames(dat3)[2] = "Gene"
colnames(dat3)[3] = "RefSeq_ID"
colnames(dat3) = gsub(".AVG_Signal", "", colnames(dat3))

## Check to see that expression data columns and sample information rows are in the same order:

all.equal(as.character(gsub("X", "", colnames(dat3)[6:101])), as.character(datSample$CoreLabel))
# [1] "60 string mismatches"

## The samples were not in the same order so reorder datSample by CoreLabel:

datSample = datSample[order(datSample[, 15]), ]

all.equal(as.character(gsub("X", "", colnames(dat3)[6:101])), as.character(datSample$CoreLabel))
# [1] TRUE

## Relabel the columns in the expression file with EC1 SampleLabel
colnames(dat3)[6:101] = as.character(datSample$SampleLabel)

## It would be nice to reorder the samples sequentially:
datSample = data.frame(datSample, Order = c(1:96))
datSample = datSample[order(datSample$Section_number), ]
datSample = data.frame(datSample, New.Order = c(1:96))
datSample = datSample[order(datSample$Order), ]

all.equal(as.character(colnames(dat3)[6:101]), as.character(datSample$SampleLabel))
# [1] TRUE

## Now reorder the column names in the expression data to match the SampleInfo.
dat3Sub = dat3[, 6:101]
dat3Sub = dat3Sub[, order(datSample$New.Order)]
dat3 = data.frame(dat3[, 1:5], dat3Sub)

datSample = datSample[order(datSample$New.Order), ]
## Check that columns are in the same order as in the sample information.
all.equal(as.character(colnames(dat3)[6:101]), as.character(datSample$SampleLabel))
# [1] TRUE

dim(dat3)
# [1] 47231   101

## Add the probe quality.

datAnnot = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1/ACC1_ExpressionData_for_SampleNetwork.csv")

datAnnot = datAnnot[order(datAnnot$PROBE_ID), ]

all.equal(datAnnot$PROBE_ID, dat3$PROBE_ID)
# [1] TRUE

dat3 = data.frame(dat3[, 1:2], ProbeQuality = datAnnot$ProbeQuality, dat3[, 3:101])

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC1")

dim(dat3)
# [1] 47231   102

## Remove the order columns from datSample
datSample = datSample[, 1:18]


## I had problems with quote being present in the data. When I write out the data table I lose some rows. To overcome this I will use quote=FALSE when I write the file out.
write.table(dat3, file = "EC1_ExpressionData_for_SampleNetwork.csv", row.names = FALSE, col.names = TRUE, sep = ",")
write.table(datSample, file = "EC1_SampleInfo_for_SampleNetwork.csv", row.names = FALSE, col.names = TRUE, sep = ",")



## Now we can run the SampleNetwork function. I did this on my local machine because it is easier with the interactivity:

rm(list = ls())

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC1")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1 = read.csv("EC1_ExpressionData_for_SampleNetwork.csv")
datSample = read.csv("EC1_SampleInfo_for_SampleNetwork.csv")

dim(dat1)
# [1] 47231   102
dim(datSample)
# [1] 96 18

indexAll = c(7:102)

SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 13,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 17,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 13, 14, 18, 11, 8, 5, 6),
    trait1 = NULL,
    asfactors1 = c(13, 14, 18),
    projectname1 = "EC1",
    cexlabels = 0.7,
    normalize1 = TRUE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)


## After the first round I removed <-2 which removed the two reference samples. No further samples removed. After quantile normalization there was a batch effect associated with ArrayID. This was corrected for with comBat using no covariates.

## After running combat there is a significant batch effect associated with ArrayPosition. To correct for this I will run another round of comBat to remove it.

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC1/EC1_SampleNetworks/ALL_02-58-50")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1 = read.csv("EC1_ALL_94_ComBat.csv")
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC1/EC1_SampleInfo_for_SampleNetwork.csv")
## Remove two reference samples from sample information file:
datSample = datSample[-c(grep("REF", datSample$SampleLabel)), ]

## Run SampleNetwork as before, except set normalize1=FALSE and subgroup1=14 to highlight array position.

indexAll = c(7:100)

SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 14,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 17,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 13, 14, 18, 11, 8, 5, 6),
    trait1 = NULL,
    asfactors1 = c(13, 14, 18),
    projectname1 = "EC1_rd2",
    cexlabels = 0.7,
    normalize1 = FALSE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)

## No samples were removed and ComBat was run on ArrayPosition with no covariates. After correcting for array position there didn't appear to be any batch effects.

## ACC2.

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC2")

## Note: added array info and position info from core spreadsheet to sample information file.  Reordered by Array_position column (e.g. 6303256015_A).

datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_SampleInfo/ACC2_SampleInfo_for_SampleNetwork.csv")

datSample = datSample[, -17]
## Add group column to datSample
datSample = data.frame(datSample, Group = "ALL")

## To translate the continuous first randomization number (QC_number) to a QC_batch:

QC_batch = ceiling((as.integer(datSample$QC_number) / 12))
datSample = data.frame(datSample, QC_batch)
## Write the amended table to the working directory.
write.table(datSample, file = "ACC2_SampleInfo_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)

## Note: open "2012-028 sample probe profile.txt" in Excel and save it as a .csv file
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2012-028 sample probe profile.csv")
dim(dat1)
# [1] 47231   398

## We will select the expression values ("AVG_Signal") columns and summarize the Detection P-values:
## We will also retain probe metadata.

seq1 = seq(3, 383, 4)

AvgDetPval = apply(dat1[, seq1 + 1], 1, mean)

## What is the distribution of mean detection P-values?

hist(AvgDetPval)

## So ~20K/47K transcripts show a mean detection P-value < .05 (i.e. are expressed significantly above background levels).

## Now we ask what fraction of transcripts are detected at P < .05 in at least one sample:

CountDetPval = apply(dat1[, seq1 + 1], 1, function(x) {
    length(x[x < .05])
})

## What is the distribution of the number of samples for which a given transcript was detected above background levels?

hist(CountDetPval, breaks = 50)

dat2 = data.frame(dat1[, 1:2], dat1[, 394], AvgDetPval, CountDetPval, dat1[, seq1])
colnames(dat2)[2] = "Gene"
colnames(dat2)[3] = "RefSeq_ID"
colnames(dat2) = gsub(".AVG_Signal", "", colnames(dat2))

## Now reorder datSample by CoreLabel:

datSample = datSample[order(datSample$CoreLabel), ]

## Give the reference samples a section number so that they can be ordered correctly:
datSample$SampleLabel = as.character(datSample$SampleLabel)

datSample[datSample$SampleLabel == "ACC2_Ref1", 16] = "REF_1"
datSample[datSample$SampleLabel == "ACC2_Ref2", 16] = "REF_2"


datSample[datSample$SampleLabel == "REF_1", "Section_number"] = 121
datSample[datSample$SampleLabel == "REF_2", "Section_number"] = 122

## Check to see that expression data columns and sample information rows are in the same order:

all.equal(as.character(gsub("X", "", colnames(dat2)[6:101])), as.character(datSample$CoreLabel))
# [1] TRUE

## So all samples are in the same order.  Lets re-label the column headers in the expression data:

colnames(dat2)[6:101] = as.character(datSample$SampleLabel)

## Now get the columns in the correct order:

datSample = datSample[order(datSample$Section_number), ]
datSample = data.frame(datSample, new.order = c(1:96))
datSample = datSample[order(datSample$CoreLabel), ]

all.equal(as.character(colnames(dat2)[6:101]), as.character(datSample$SampleLabel))
# [1] TRUE

dat2Sub = dat2[, 6:101]
dat2Sub = dat2Sub[, order(datSample$new.order)]
datSample = datSample[order(datSample$Section_number), ]

all.equal(as.character(colnames(dat2Sub)), as.character(datSample$SampleLabel))
# [1] TRUE

dat2 = data.frame(dat2[, 1:5], dat2Sub)

## Reorder the expression data by PROBE_ID:

dat2 = dat2[order(dat2$PROBE_ID), ]

## Add the probe quality.

datAnnot = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1/ACC1_ExpressionData_for_SampleNetwork.csv")

datAnnot = datAnnot[order(datAnnot$PROBE_ID), ]

all.equal(datAnnot$PROBE_ID, dat2$PROBE_ID)
# [1] TRUE

dat2 = data.frame(dat2[, 1:2], ProbeQuality = datAnnot$ProbeQuality, dat2[, 3:101])

write.table(dat2, file = "ACC2_ExpressionData_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)
write.table(datSample, file = "ACC2_SampleInfo_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)

rm(list = ls())

## Run SampleNetwork:

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC2")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1 = read.csv("ACC2_ExpressionData_for_SampleNetwork.csv")
datSample = read.csv("ACC2_SampleInfo_for_SampleNetwork.csv")

indexAll = c(7:102)

## Note: coloring by ArrayID to start (subgroup1=13), which may or may not be a significant batch effect.

SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 13,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 17,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 13, 14, 18, 11, 8, 5, 6),
    trait1 = NULL,
    asfactors1 = c(13, 14, 18),
    projectname1 = "ACC2",
    cexlabels = 0.7,
    normalize1 = TRUE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)


## After the first round I removed samples <-2. Second round I removed samples <-3 and after the third round <-3 and <-3 after the 4th round. The fifth round I removed <-3 which removed 1 sample. A total of 13 samples were removed. After quantile normalization there was a batch effect associated with ArrayID which I corrected for with comBat using no covariates.



## EC2


setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC2")

## Note: added array info and position info from core spreadsheet to sample information file.  Reordered by Array_position column (e.g. 6303256015_A).

datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_SampleInfo/EC2_SampleInfo_for_SampleNetwork.csv")

datSample = datSample[, -c(17:18)]
## Add group column to datSample
datSample = data.frame(datSample, Group = "ALL")

## To translate the continuous first randomization number (QC_number) to a QC_batch:

QC_batch = ceiling((as.integer(datSample$QC_number) / 12))
datSample = data.frame(datSample, QC_batch)
## Write the amended table to the working directory.
write.table(datSample, file = "EC2_SampleInfo_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)

## Note: open "2012-125 sample probe profile.txt" in Excel and save it as a .csv file
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2012-125 sample probe profile.csv")
dim(dat1)
# [1] 47231   398

## We will select the expression values ("AVG_Signal") columns and summarize the Detection P-values:
## We will also retain probe metadata.

seq1 = seq(3, 383, 4)

AvgDetPval = apply(dat1[, seq1 + 1], 1, mean)

## What is the distribution of mean detection P-values?

hist(AvgDetPval)

## So ~20K/47K transcripts show a mean detection P-value < .05 (i.e. are expressed significantly above background levels).

## Now we ask what fraction of transcripts are detected at P < .05 in at least one sample:

CountDetPval = apply(dat1[, seq1 + 1], 1, function(x) {
    length(x[x < .05])
})

## What is the distribution of the number of samples for which a given transcript was detected above background levels?

hist(CountDetPval, breaks = 50)

dat2 = data.frame(dat1[, 1:2], dat1[, 394], AvgDetPval, CountDetPval, dat1[, seq1])
colnames(dat2)[2] = "Gene"
colnames(dat2)[3] = "RefSeq_ID"
colnames(dat2) = gsub(".AVG_Signal", "", colnames(dat2))

## Now reorder datSample by CoreLabel:

datSample = datSample[order(datSample$CoreLabel), ]

## Give the reference samples a section number so that they can be ordered correctly:
datSample$SampleLabel = as.character(datSample$SampleLabel)

datSample[datSample$SampleLabel == "EC2_Ref1", 16] = "REF_1"
datSample[datSample$SampleLabel == "EC2_Ref2", 16] = "REF_2"


datSample[datSample$SampleLabel == "REF_1", "Section_number"] = 132
datSample[datSample$SampleLabel == "REF_2", "Section_number"] = 133

## Check to see that expression data columns and sample information rows are in the same order:

all.equal(as.character(gsub("X", "", colnames(dat2)[6:101])), as.character(datSample$CoreLabel))
# [1] TRUE

## So all samples are in the same order.  Lets re-label the column headers in the expression data:

colnames(dat2)[6:101] = as.character(datSample$SampleLabel)

## Now get the columns in the correct order:

datSample = datSample[order(datSample$Section_number), ]
datSample = data.frame(datSample, new.order = c(1:96))
datSample = datSample[order(datSample$CoreLabel), ]

all.equal(as.character(colnames(dat2)[6:101]), as.character(datSample$SampleLabel))
# [1] TRUE

dat2Sub = dat2[, 6:101]
dat2Sub = dat2Sub[, order(datSample$new.order)]
datSample = datSample[order(datSample$Section_number), ]

all.equal(as.character(colnames(dat2Sub)), as.character(datSample$SampleLabel))
# [1] TRUE

dat2 = data.frame(dat2[, 1:5], dat2Sub)

## Reorder the expression data by PROBE_ID:

dat2 = dat2[order(dat2$PROBE_ID), ]

## Add the probe quality.

datAnnot = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1/ACC1_ExpressionData_for_SampleNetwork.csv")

datAnnot = datAnnot[order(datAnnot$PROBE_ID), ]

all.equal(datAnnot$PROBE_ID, dat2$PROBE_ID)
# [1] TRUE

dat2 = data.frame(dat2[, 1:2], ProbeQuality = datAnnot$ProbeQuality, dat2[, 3:101])

write.table(dat2, file = "EC2_ExpressionData_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)
write.table(datSample, file = "EC2_SampleInfo_for_SampleNetwork.csv", sep = ",", row.names = F, col.names = T)

rm(list = ls())

## Run SampleNetwork (SampleNetwork run locally because it is easier to run interactively than on the cluster):

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC2")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1 = read.csv("EC2_ExpressionData_for_SampleNetwork.csv")
datSample = read.csv("EC2_SampleInfo_for_SampleNetwork.csv")

indexAll = c(7:102)

## Note: coloring by ArrayID to start (subgroup1=13), which may or may not be a significant batch effect.

SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 13,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 17,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 13, 14, 18, 11, 8, 5, 6),
    trait1 = NULL,
    asfactors1 = c(13, 14, 18),
    projectname1 = "EC2",
    cexlabels = 0.7,
    normalize1 = TRUE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)

## After round 1 I chose <-2 which removed the two reference samples. After round 2 I removed <-3 which removed another 2 samples. A total of 4 samples were removed, 2 reference samples and 2 actual samples. After quantile normalization there was a significant batch effects associated with QC_Batch and ArrayID. I corrected QC_Batch using comBat with no covariates.

## After using comBat to correct for QC_Batch ArrayID was still a significant batch effect so I ran a second round of comBat to correct for ArrayID. Normalize was set to FALSE

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC2/EC2_SampleNetworks/ALL_03-20-14")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC2/EC2_SampleNetworks/ALL_03-20-14/EC2_ALL_92_ComBat.csv")
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC2/EC2_SampleInfo_for_SampleNetwork.csv")

datout = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC2/EC2_SampleNetworks/EC2_Group_outliers_03-20-14.csv")


## Remove outliers from the sampleinfo file:
datSample = datSample[!is.element(datSample$SampleLabel, datout$SampleLabel), ]

## Run SampleNetwork as before, except set normalize1=FALSE and subgroup1=14 to highlight array position.

indexAll = c(7:98)

SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 14,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 17,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 13, 14, 18, 11, 8, 5, 6),
    trait1 = NULL,
    asfactors1 = c(13, 14, 18),
    projectname1 = "EC2_rd2",
    cexlabels = 0.7,
    normalize1 = FALSE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)

## No further samples were removed and comBat was run to correct for ArrayID with no covariates.
