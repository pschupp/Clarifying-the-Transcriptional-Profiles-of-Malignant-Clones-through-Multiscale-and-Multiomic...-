setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor")
## Read in the sample information file:
datSample = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_SampleInfo/SF9495_SampleInfo_for_SampleNetwork.csv")
colnames(datSample)
# [1] "Original"            "Randomized1"         "Randomized2"
# [4] "ND_RNA_conc_ng_ul"   "ND_RNA_260.280"      "ND_RNA_260.230"
# [7] "BA__RNA_conc_ng_ul"  "RIN"                 "Ratio_28s.18s"
# [10] "ND_DNA_conc_ng_ul"   "ND_DNA_260.280"      "ND_DNA_260.230"
# [13] "Mean.RNA.Conc.ng_ul" "FirstScreen"
## I will rename some of the columns:
colnames(datSample)[1] = "section_number"
colnames(datSample)[2] = "QC_number"
colnames(datSample)[3] = "Core_number"
## Add a samplelabel:
SampleLabel = paste("SF9495_", datSample$section_number, sep = "")
datSample = data.frame(datSample, SampleLabel = SampleLabel)
dim(datSample)
# [1] 72 15
## Create some custom names for the reference samples and some section numbers so that I can sort them:
datSample[, 15] = as.character(datSample[, 15])
datSample[71, 15] = "SF9495_REF_1"
datSample[72, 15] = "SF9495_REF_2"
datSample[71, 1] = 82
datSample[72, 1] = 83
## It would also be useful to add the array ID and position so I can correct for any batch effects. I opened the sample key file that I got back from the UCLA core and saved it as a csv file:
datkey = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2013-140_Sample Key.csv")
## Rename the reference samples within the sample_key information file:
datkey[, 4] = as.character(datkey[, 4])
datkey[1, 4] = "REF_1"
datkey[72, 4] = "REF_2"
datkey[, 4] = paste("SF9495", datkey[, 4], sep = "_")
## These samples should be in the same order as in the expression file. I will make an order vector so that I can resort these and the sample information file.
datkey = data.frame(datkey, Order = c(1:72))
datkey = datkey[order(datkey[, 4]), ]
datSample = datSample[order(datSample$SampleLabel), ]
## The sampleinfo and samplekey form UCLA should be in the same order now.
## Check that the samples are in the correct order:
all.equal(as.character(datSample[, 15]), as.character(datkey[, 4]))
# [1] TRUE
datSample = data.frame(datSample, Array_ID = datkey$general.array, Array_Position = datkey$genexstripe.controling.stripe, group = "ALL")
QC_Batch = ceiling(datSample$QC_number / 12)
datSample = data.frame(datSample[, 1:2], QC_Batch = QC_Batch, datSample[, 3:18])
## Now I will load the expression file and relabel the column names. I opened the file 2013-140_Sample probe profile.xlsx and resaved it as csv:
dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Raw_microarray_data/2013-140 sample probe profile.csv")
dim(dat1)
# [1] 47202   302
## Take the average signal for each of these samples:
seq1 = seq(3, 287, 4)
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
dat2 = data.frame(dat1[, 1:2], dat1[, 298], AvgDetPval, CountDetPval, dat1[, seq1])
colnames(dat2)[2] = "Gene"
colnames(dat2)[3] = "RefSeq_ID"
colnames(dat2) = gsub(".AVG_Signal", "", colnames(dat2))
## Now I will create a Core_ID by pasting together the Array_ID and Array_Position. I can use this to rename the columns in the expression data:
Core_ID = paste("X", datSample$Array_ID, "_", datSample$Array_Position, sep = "")
datSample = data.frame(datSample, Core_ID = Core_ID)
## Check that the samples are in the same order between the sample information and expression files:
datSample = datSample[order(datSample$Core_ID), ]
all.equal(as.character(datSample$Core_ID), as.character(colnames(dat2)[6:77]))
# [1] TRUE
## Now I will rename the columns in dat2 with datSample$SampleLabel:
colnames(dat2)[6:77] = as.character(datSample$SampleLabel)
## Now I will reorder the sample information file and the expression file by section number:
datSample = datSample[order(datSample$section_number), ]
datSample = data.frame(datSample, new.order = c(1:length(datSample[, 1])))
datSample = datSample[order(datSample$Core_ID), ]
all.equal(as.character(datSample$SampleLabel), as.character(colnames(dat2)[6:77]))
# [1] TRUE
## Now reorder the expression file:
datSub = dat2[, 6:77]
all.equal(as.character(datSample$SampleLabel), as.character(colnames(datSub)))
# [1] TRUE
datSub = datSub[, order(datSample$new.order)]
dat2 = data.frame(dat2[, 1:5], datSub)
datSample = datSample[order(datSample$section_number), ]
## I also want to get the probe annotations so I will read in the data frame from the ACC1 and copy them across:
dat3 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1/ACC1_ExpressionData_for_SampleNetwork.csv")
dat2 = dat2[order(dat2$PROBE_ID), ]
dat3 = dat3[order(dat3$PROBE_ID), ]
dim(dat2)
# [1] 47202    77
dim(dat3)
# [1] 47231    98
## There are more probes for ACC1 because some probes were excluded by the UCLA core for the glioma SF9495 dataset. I will take just the intersect with the glioma dataset:
dat3 = dat3[is.element(dat3$PROBE_ID, intersect(dat3$PROBE_ID, dat2$PROBE_ID)), ]
dim(dat3)
# [1] 47202    98
all.equal(as.character(dat2$PROBE_ID), as.character(dat3$PROBE_ID))
# [1] TRUE
dat2 = data.frame(dat2[, 1:5], ProbeQuality = dat3$ProbeQuality, dat2[6:77])
## Now write the sample information and expression files out:
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor")
write.table(datSample, file = "SF9495_SampleInfo_for_SampleNetwork.csv", col.names = TRUE, row.names = FALSE, sep = ",")
write.table(dat2, file = "SF9495_ExpressionFile_for_SampleNetwork.csv", col.names = TRUE, row.names = FALSE, sep = ",")
## Now to run the SampleNetwork Function. I ran this locally due to interactivity needed:
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor")
dat1 = read.csv("SF9495_ExpressionFile_for_SampleNetwork.csv")
datSample = read.csv("SF9495_SampleInfo_for_SampleNetwork.csv")
source("/fast-data/Shared/Code/SJS/SampleNetwork_1.01_cluster_SJS.r")
indexAll = c(7:78)
## Color by QC batch (col3):
SampleNetwork(
    datExprT = dat1,
    method1 = "correlation",
    impute1 = FALSE,
    subset1 = NULL,
    skip1 = 6,
    indices1 = list(indexAll),
    subgroup1 = 3,
    sampleinfo1 = datSample,
    samplelabels1 = 16,
    grouplabels1 = 19,
    fitmodels1 = TRUE,
    whichmodel1 = "univariate",
    whichfit1 = "pc1",
    btrait1 = c(1, 3, 5, 6, 7, 8, 9, 10, 14, 17, 18),
    trait1 = NULL,
    asfactors1 = c(3, 9, 17, 18),
    projectname1 = "SF9495",
    cexlabels = 0.7,
    normalize1 = TRUE,
    replacenegs1 = FALSE,
    exportfigures1 = TRUE,
    verbose = TRUE
)
## After the first round the data looked pretty good and the two reference samples looked clearly different to the rest of the samples. I removed these by using <-2.
## After the second round there was one outlier (section 48) which I removed using <-3:
## No further samples removed, 3 samples in total were removed, the controls and 1 other sample:
## There were no further significant batch effects other than section number but since this represents variation between samples this is what I would expect and it will not be corrected for.
## After removing 1 outlier which was sample 48 this leaves me with 69 samples.
