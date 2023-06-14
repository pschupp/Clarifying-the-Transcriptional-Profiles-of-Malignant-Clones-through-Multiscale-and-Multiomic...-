## use the normal brain data as a reference to determine the copy number.
# {{{
.libPaths("/fast-data/R-3.1.1-packages-SJS")
library(ChAMP)
## Set the directory to the location of the idat files:
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Methylation_analysis_033116/2015-9159_iScans_Combined")
## Open my R data object with raw intensity values for my data
load("myLoad.RData")
## Now load the brain control samples:
Ctls = read.csv("/big-data2/Sam/GSE43414_Brain_450K_SJS/GSE43414_signal_intensity_controls_only.csv")
myLoad$pd = data.frame(Sample_Name = myLoad$pd$Sample_Name, Sample_Plate = myLoad$pd$Sample_Plate, Sample_Group = myLoad$pd$Sample_Group, Pool_ID = myLoad$pd$Pool_ID, Project = NA, Sample_Well = myLoad$pd$Sample_Well, Array = myLoad$pd$Array, Slide = myLoad$pd$Slide, Basename = myLoad$pd$Basename, filenames = myLoad$pd$filenames)
myLoad$pd$Sample_Group = "SF9495"
## I now need to get the sampleinformation for the control samples in the same format.
datSample = read.csv("/big-data2/Sam/GSE43414_Brain_450K_SJS/GSE43414_sampleinfo_controls_only.csv")
datSample = data.frame(Sample_Name = datSample$barcode, Sample_Plate = "Plate2", Sample_Group = "Control", Pool_ID = NA, Project = NA, Sample_Well = NA, Array = datSample$array.pos, Slide = datSample$array.id, Basename = NA, filenames = NA)
## Check that the samples are in the correct order between the sample info and intensity data:
all.equal(as.character(datSample$Sample_Name), as.character(colnames(Ctls)[2:95]))
# [1] TRUE
all.equal(as.character(myLoad$pd$Sample_Name), as.character(colnames(myLoad$intensity)))
# [1] TRUE
## Get the arrays in the same order:
myLoad$intensity = myLoad$intensity[order(rownames(myLoad$intensity)), ]
rownames(Ctls) = Ctls[, 1]
Ctls = Ctls[, 2:95]
Ctls = Ctls[order(rownames(Ctls)), ]
dim(Ctls)
# [1] 485577     94
dim(myLoad$intensity)
# [1] 461797     69
## Subset the ctrls down to the same probes as in the expression data:
Ctls = Ctls[is.element(rownames(Ctls), rownames(myLoad$intensity)), ]
dim(Ctls)
# [1] 461797     94
all.equal(as.character(rownames(Ctls)), as.character(rownames(myLoad$intensity)))
# [1] TRUE
## I will run the code from champ.CNA manually so that I can get a logR ratio.
## Set pd for my samples:
pd = myLoad$pd

data(probe.features)
normalize.quantiles <- NULL
rm(normalize.quantiles)
control.intsqnlog <- NULL
CNA <- NA
rm(CNA)
smooth.CNA <- NA
rm(smooth.CNA)
segment <- NA
rm(segment)

## Set the ints as intensity file:
ints = myLoad$intensity
## Now set the ctrl intensities:
ctlIntensity = Ctls
ctlIntensity = ctlIntensity[which(row.names(ctlIntensity) %in% row.names(ints)), ]
all.equal(rownames(ctlIntensity), rownames(ints))
# [1] TRUE
ints = cbind(ints, ctlIntensity)
all.equal(colnames(pd), colnames(datSample))
# [1] TRUE
pd = rbind(pd, datSample)
## Check that everything is in the same order:
all.equal(as.character(pd$Sample_Name), colnames(ints))
# [1] TRUE
## I need preprocessCore for the normalize.quantiles function
.libPaths("/fast-data/R-3.1.1-packages-SJS")
library(preprocessCore)
library(DNAcopy)
names <- colnames(ints)
## The combined normal
intsqn <- normalize.quantiles(as.matrix(ints))
colnames(intsqn) <- names
intsqnlog <- log2(intsqn)

## Now define the controlSamples and caseSamples
controlGroup = "Control"
controlSamples = pd[which(pd$Sample_Group == controlGroup), ]
caseSamples = pd[which(pd$Sample_Group != controlGroup), ]
##  Now I want to log2 transform the intensity values:
case.intsqnlog <- intsqnlog[, which(colnames(intsqnlog) %in% caseSamples$Sample_Name)]
control.intsqnlog <- intsqnlog[, which(colnames(intsqnlog) %in% controlSamples$Sample_Name)]
## Calculate the mean log2 intensity values for the controls
control.intsqnlog <- rowMeans(control.intsqnlog)
## Create a new table with the log intensity values for the tumor samples:
intsqnlogratio <- case.intsqnlog
## Now calculate the logR ratios:
for (i in 1:ncol(case.intsqnlog)) {
    intsqnlogratio[, i] <- case.intsqnlog[, i] - control.intsqnlog
}
## Now I would like to add back the probe information:
library(ChAMP)
data(probe.features)
ints <- data.frame(ints, probe.features$MAPINFO[match(rownames(ints), rownames(probe.features))])
names(ints)[length(ints)] <- "MAPINFO"
ints <- data.frame(ints, probe.features$CHR[match(rownames(ints), rownames(probe.features))])
names(ints)[length(ints)] <- "CHR"
levels(ints$CHR)[levels(ints$CHR) == "X"] = "23"
levels(ints$CHR)[levels(ints$CHR) == "Y"] = "24"
CHR <- as.numeric(levels(ints$CHR))[ints$CHR]
ints$MAPINFO <- as.numeric(ints$MAPINFO)
MAPINFO = probe.features$MAPINFO[match(rownames(ints), rownames(probe.features))]
MAPINFO <- as.numeric(MAPINFO)
## At this point I could create a table of logR ratios:
logRdf = data.frame(chrom = CHR, pos = MAPINFO, intsqnlogratio)
logRdf = logRdf[order(logRdf$chrom, logRdf$pos), ]
dim(logRdf)
# [1] 461797     71
## Write this file out for future reference:
write.table(logRdf, file = "LogR_ratio_tumor_vs_normal_brain.csv", col.names = TRUE, row.names = FALSE, sep = ",")
library(copynumber)
## Also try running multipcf with normalize set to FALSE since the intensity values have already been quantile normalized.
multi.seg <- multipcf(data = logRdf, verbose = FALSE, normalize = FALSE, assembly = "hg19", gamma =20)
dim(multi.seg)
# [1] 40 74
## Order the segmentation file by chromosome and then position
multi.seg$chrom = as.numeric(multi.seg$chrom)
multi.seg = multi.seg[order(multi.seg$chrom, multi.seg$start.pos), ]
## There are far fewer copy number changes when using the normal brain to normalize to rather than blood. Potentially the additional CNVs that I saw when comparing to blood reflect differences in methylation between brain and blood rather than differences between tumor and normal brain.
write.table(multi.seg, file = "Segmentation_of_all_sections_copynumberR_relative_to_normal_brain.csv", col.names = TRUE, row.names = FALSE, sep = ",")
# }}} 
