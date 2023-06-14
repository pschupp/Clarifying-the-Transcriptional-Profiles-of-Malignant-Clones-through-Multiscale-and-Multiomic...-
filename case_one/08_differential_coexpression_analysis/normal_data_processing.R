## Build subtraction matrices using my new GBM dataset and subtracting the other ACC and EC datasets. Summarize the 4 different subtraction matrices by taking the parallel minimum for all of the 4 subtractions.

## First I will build the individual subtraction matrices:

.libPaths("/fast-data/R-3.1.1-packages")
library(WGCNA)

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")
dat2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC1/Renumbered_ACC1_ALL_92_ComBat.csv")

dat1 = dat1[order(dat1$PROBE_ID), ]
dat2 = dat2[order(dat2$PROBE_ID), ]

dim(dat1)
# [1] 47202    75
dim(dat2)
# [1] 47231    98

## They are different numbers of probes:

dat2 = dat2[is.element(dat2$PROBE_ID, intersect(dat2$PROBE_ID, dat1$PROBE_ID)), ]

dim(dat2)
# [1] 47202    98

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

## I need to write out the subset file and reread it in otherwise R still sees it at 47231 in CorMatDNA.
write.table(dat2, file = "ACC1_Expression_ProbesMatchedto_SF9495.csv", col.names = TRUE, row.names = FALSE, sep = ",")
dat2 = read.csv("ACC1_Expression_ProbesMatchedto_SF9495.csv")

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

dim(dat2)
# [1] 47202    98


## Read in the new CorMatDNA function from Mike:

CorMatDNA = function(dat1, sampleindex1, simType, name1, WD1, dat2, sampleindex2, name2, WD2, subset, difference = c("dat1-dat2", "dat2-dat1", "reciprocal")) {
    timestamp()

    rtoz = function(x) {
        0.5 * log((1 + x) / (1 - x))
    }

    ztor = function(x) {
        (exp(2 * x) - 1) / (exp(2 * x) + 1)
    }

    if (class(dat1[, 1]) == "integer" | class(dat1[, 1]) == "numeric") {
        dat1 = dat1[order(as.numeric(as.character(dat1[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat1[, 1])))]
    } else {
        dat1 = dat1[order(as.character(dat1[, 1])), ]
        subset = subset[order(as.character(dat1[, 1]))]
    }

    if (class(dat2[, 1]) == "integer" | class(dat2[, 1]) == "numeric") {
        dat2 = dat2[order(as.numeric(as.character(dat2[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat2[, 1])))]
    } else {
        dat2 = dat2[order(as.character(dat2[, 1])), ]
        subset = subset[order(as.character(dat2[, 1]))]
    }

    if (all(dat1[, 1] == dat2[, 1]) == FALSE) {
        stop("dat1[,1] and dat2[,1] are not equivalent")
    }

    if (!is.null(subset)) {
        subset = as.logical(subset)
        datExpr1 = t(dat1[subset, sampleindex1])
        colnames(datExpr1) = as.character(dat1[subset, 1])
        datExpr2 = t(dat2[subset, sampleindex2])
        colnames(datExpr2) = as.character(dat2[subset, 1])
    } else {
        datExpr1 = t(dat1[, sampleindex1])
        colnames(datExpr1) = as.character(dat1[, 1])
        datExpr2 = t(dat2[, sampleindex2])
        colnames(datExpr2) = as.character(dat2[, 1])
    }

    if (simType == "Pearson") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "p", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "p", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Spearman") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "s", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "s", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Bicor") {
        print("Calculating corMat1...")
        corMat1 = bicor(datExpr1, use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = bicor(datExpr2, use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    # print("Normalizing corMat1...")
    # corMat1=((ztor(scale(rtoz(corMat1))))+1)/2
    # corMat1=(ztor(rtoz(corMat1)/(1/sqrt(length(sampleindex1)-3)))+1)/2
    # diag(corMat1)=0
    # collectGarbage()

    if (min(corMat1, na.rm = T) < (-1) | max(corMat1, na.rm = T) > 1) {
        stop("Correlations in corMat1 outside of {-1,1}")
    }

    # print("Normalizing corMat2...")
    # corMat2=((ztor(scale(rtoz(corMat2))))+1)/2
    # corMat2=(ztor(rtoz(corMat2)/(1/sqrt(length(sampleindex2)-3)))+1)/2
    # diag(corMat2)=0
    # collectGarbage()

    if (min(corMat2, na.rm = T) < (-1) | max(corMat2, na.rm = T) > 1) {
        stop("Correlations in corMat2 outside of {-1,1}")
    }

    if (difference == "reciprocal") {
        simMat = corMat2 - corMat1
        filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        setwd(WD2)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
        simMat = corMat1 - corMat2
        filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        setwd(WD1)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
    } else {
        if (difference == "dat2-dat1") {
            setwd(WD2)
            simMat = corMat2 - corMat1
            filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        }

        if (difference == "dat1-dat2") {
            setwd(WD1)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")

            setwd(WD2)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        }

        collectGarbage()
        save(simMat, file = filenameout)
    }

    collectGarbage()

    timestamp()
} ## end of function


CorMatDNA(
    dat1 = dat1,
    sampleindex1 = c(7:75),
    name1 = "SF9495_PandG_30k",
    WD1 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    dat2 = dat2,
    sampleindex2 = c(7:98),
    name2 = "ACC1_PandG_30k",
    WD2 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    subset = !is.na(dat1$ProbeQuality) & (dat1$ProbeQuality == "P" | dat1$ProbeQuality == "G"),
    difference = "reciprocal",
    simType = "Bicor"
)


## Now do the same for ACC2:

.libPaths("/fast-data/R-3.1.1-packages")
library(WGCNA)

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")
dat2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/ACC2/Renumbered_ACC2_ALL_83_ComBat.csv")

dat1 = dat1[order(dat1$PROBE_ID), ]
dat2 = dat2[order(dat2$PROBE_ID), ]

dim(dat1)
# [1] 47202    75
dim(dat2)
# [1] 47231    89

## They are different numbers of probes:

dat2 = dat2[is.element(dat2$PROBE_ID, intersect(dat2$PROBE_ID, dat1$PROBE_ID)), ]

dim(dat2)
# [1] 47202    89

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

## Write out the subset table and reread it into R so that corMatDNA sees it as the correct dimensions.
write.table(dat2, file = "ACC2_Expression_ProbesMatchedto_SF9495.csv", col.names = TRUE, row.names = FALSE, sep = ",")

dat2 = read.csv("ACC2_Expression_ProbesMatchedto_SF9495.csv")

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

dim(dat2)
# [1] 47202    89

## Read in the new CorMatDNA function from Mike:

CorMatDNA = function(dat1, sampleindex1, simType, name1, WD1, dat2, sampleindex2, name2, WD2, subset, difference = c("dat1-dat2", "dat2-dat1", "reciprocal")) {
    timestamp()

    rtoz = function(x) {
        0.5 * log((1 + x) / (1 - x))
    }

    ztor = function(x) {
        (exp(2 * x) - 1) / (exp(2 * x) + 1)
    }

    if (class(dat1[, 1]) == "integer" | class(dat1[, 1]) == "numeric") {
        dat1 = dat1[order(as.numeric(as.character(dat1[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat1[, 1])))]
    } else {
        dat1 = dat1[order(as.character(dat1[, 1])), ]
        subset = subset[order(as.character(dat1[, 1]))]
    }

    if (class(dat2[, 1]) == "integer" | class(dat2[, 1]) == "numeric") {
        dat2 = dat2[order(as.numeric(as.character(dat2[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat2[, 1])))]
    } else {
        dat2 = dat2[order(as.character(dat2[, 1])), ]
        subset = subset[order(as.character(dat2[, 1]))]
    }

    if (all(dat1[, 1] == dat2[, 1]) == FALSE) {
        stop("dat1[,1] and dat2[,1] are not equivalent")
    }

    if (!is.null(subset)) {
        subset = as.logical(subset)
        datExpr1 = t(dat1[subset, sampleindex1])
        colnames(datExpr1) = as.character(dat1[subset, 1])
        datExpr2 = t(dat2[subset, sampleindex2])
        colnames(datExpr2) = as.character(dat2[subset, 1])
    } else {
        datExpr1 = t(dat1[, sampleindex1])
        colnames(datExpr1) = as.character(dat1[, 1])
        datExpr2 = t(dat2[, sampleindex2])
        colnames(datExpr2) = as.character(dat2[, 1])
    }

    if (simType == "Pearson") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "p", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "p", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Spearman") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "s", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "s", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Bicor") {
        print("Calculating corMat1...")
        corMat1 = bicor(datExpr1, use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = bicor(datExpr2, use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    # print("Normalizing corMat1...")
    # corMat1=((ztor(scale(rtoz(corMat1))))+1)/2
    # corMat1=(ztor(rtoz(corMat1)/(1/sqrt(length(sampleindex1)-3)))+1)/2
    # diag(corMat1)=0
    # collectGarbage()

    if (min(corMat1, na.rm = T) < (-1) | max(corMat1, na.rm = T) > 1) {
        stop("Correlations in corMat1 outside of {-1,1}")
    }

    # print("Normalizing corMat2...")
    # corMat2=((ztor(scale(rtoz(corMat2))))+1)/2
    # corMat2=(ztor(rtoz(corMat2)/(1/sqrt(length(sampleindex2)-3)))+1)/2
    # diag(corMat2)=0
    # collectGarbage()

    if (min(corMat2, na.rm = T) < (-1) | max(corMat2, na.rm = T) > 1) {
        stop("Correlations in corMat2 outside of {-1,1}")
    }

    if (difference == "reciprocal") {
        simMat = corMat2 - corMat1
        filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        setwd(WD2)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
        simMat = corMat1 - corMat2
        filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        setwd(WD1)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
    } else {
        if (difference == "dat2-dat1") {
            setwd(WD2)
            simMat = corMat2 - corMat1
            filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        }

        if (difference == "dat1-dat2") {
            setwd(WD1)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")

            setwd(WD2)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        }

        collectGarbage()
        save(simMat, file = filenameout)
    }

    collectGarbage()

    timestamp()
} ## end of function


CorMatDNA(
    dat1 = dat1,
    sampleindex1 = c(7:75),
    name1 = "SF9495_PandG_30k",
    WD1 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    dat2 = dat2,
    sampleindex2 = c(7:89),
    name2 = "ACC2_PandG_30k",
    WD2 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    subset = !is.na(dat1$ProbeQuality) & (dat1$ProbeQuality == "P" | dat1$ProbeQuality == "G"),
    difference = "reciprocal",
    simType = "Bicor"
)




## Now I will do the same for EC1:
.libPaths("/fast-data/R-3.1.1-packages")
library(WGCNA)

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")
dat2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC1/Renumbered_EC1_rd2_ALL_94_ComBat.csv")

dat1 = dat1[order(dat1$PROBE_ID), ]
dat2 = dat2[order(dat2$PROBE_ID), ]

dim(dat1)
# [1] 47202    75
dim(dat2)
# [1] 47231   100

## They are different numbers of probes:

dat2 = dat2[is.element(dat2$PROBE_ID, intersect(dat2$PROBE_ID, dat1$PROBE_ID)), ]

dim(dat2)
# [1] 47202   100

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

write.table(dat2, file = "EC1_Expression_ProbesMatchedto_SF9495.csv", col.names = TRUE, row.names = FALSE, sep = ",")

dat2 = read.csv("EC1_Expression_ProbesMatchedto_SF9495.csv")

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

dim(dat2)
# [1] 47202   100

## Read in the new CorMatDNA function from Mike:

CorMatDNA = function(dat1, sampleindex1, simType, name1, WD1, dat2, sampleindex2, name2, WD2, subset, difference = c("dat1-dat2", "dat2-dat1", "reciprocal")) {
    timestamp()

    rtoz = function(x) {
        0.5 * log((1 + x) / (1 - x))
    }

    ztor = function(x) {
        (exp(2 * x) - 1) / (exp(2 * x) + 1)
    }

    if (class(dat1[, 1]) == "integer" | class(dat1[, 1]) == "numeric") {
        dat1 = dat1[order(as.numeric(as.character(dat1[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat1[, 1])))]
    } else {
        dat1 = dat1[order(as.character(dat1[, 1])), ]
        subset = subset[order(as.character(dat1[, 1]))]
    }

    if (class(dat2[, 1]) == "integer" | class(dat2[, 1]) == "numeric") {
        dat2 = dat2[order(as.numeric(as.character(dat2[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat2[, 1])))]
    } else {
        dat2 = dat2[order(as.character(dat2[, 1])), ]
        subset = subset[order(as.character(dat2[, 1]))]
    }

    if (all(dat1[, 1] == dat2[, 1]) == FALSE) {
        stop("dat1[,1] and dat2[,1] are not equivalent")
    }

    if (!is.null(subset)) {
        subset = as.logical(subset)
        datExpr1 = t(dat1[subset, sampleindex1])
        colnames(datExpr1) = as.character(dat1[subset, 1])
        datExpr2 = t(dat2[subset, sampleindex2])
        colnames(datExpr2) = as.character(dat2[subset, 1])
    } else {
        datExpr1 = t(dat1[, sampleindex1])
        colnames(datExpr1) = as.character(dat1[, 1])
        datExpr2 = t(dat2[, sampleindex2])
        colnames(datExpr2) = as.character(dat2[, 1])
    }

    if (simType == "Pearson") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "p", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "p", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Spearman") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "s", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "s", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Bicor") {
        print("Calculating corMat1...")
        corMat1 = bicor(datExpr1, use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = bicor(datExpr2, use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    # print("Normalizing corMat1...")
    # corMat1=((ztor(scale(rtoz(corMat1))))+1)/2
    # corMat1=(ztor(rtoz(corMat1)/(1/sqrt(length(sampleindex1)-3)))+1)/2
    # diag(corMat1)=0
    # collectGarbage()

    if (min(corMat1, na.rm = T) < (-1) | max(corMat1, na.rm = T) > 1) {
        stop("Correlations in corMat1 outside of {-1,1}")
    }

    # print("Normalizing corMat2...")
    # corMat2=((ztor(scale(rtoz(corMat2))))+1)/2
    # corMat2=(ztor(rtoz(corMat2)/(1/sqrt(length(sampleindex2)-3)))+1)/2
    # diag(corMat2)=0
    # collectGarbage()

    if (min(corMat2, na.rm = T) < (-1) | max(corMat2, na.rm = T) > 1) {
        stop("Correlations in corMat2 outside of {-1,1}")
    }

    if (difference == "reciprocal") {
        simMat = corMat2 - corMat1
        filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        setwd(WD2)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
        simMat = corMat1 - corMat2
        filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        setwd(WD1)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
    } else {
        if (difference == "dat2-dat1") {
            setwd(WD2)
            simMat = corMat2 - corMat1
            filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        }

        if (difference == "dat1-dat2") {
            setwd(WD1)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")

            setwd(WD2)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        }

        collectGarbage()
        save(simMat, file = filenameout)
    }

    collectGarbage()

    timestamp()
} ## end of function


CorMatDNA(
    dat1 = dat1,
    sampleindex1 = c(7:75),
    name1 = "SF9495_PandG_30k",
    WD1 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    dat2 = dat2,
    sampleindex2 = c(7:100),
    name2 = "EC1_PandG_30k",
    WD2 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    subset = !is.na(dat1$ProbeQuality) & (dat1$ProbeQuality == "P" | dat1$ProbeQuality == "G"),
    difference = "reciprocal",
    simType = "Bicor"
)


## Finally prepare the EC2 subtraction matrix:

.libPaths("/fast-data/R-3.1.1-packages")
library(WGCNA)

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")
dat2 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/EC2/Renumbered_EC2_rd2_ALL_92_ComBat.csv")

dat1 = dat1[order(dat1$PROBE_ID), ]
dat2 = dat2[order(dat2$PROBE_ID), ]

dim(dat1)
# [1] 47202    75
dim(dat2)
# [1] 47231    98

## They are different numbers of probes:

dat2 = dat2[is.element(dat2$PROBE_ID, intersect(dat2$PROBE_ID, dat1$PROBE_ID)), ]

dim(dat2)
# [1] 47202    98

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

## Write out the subset data frame and reread it into R so that CorMatDNA sees it as the correct dimensions
write.table(dat2, file = "EC2_Expression_ProbesMatchedto_SF9495.csv", col.names = TRUE, row.names = FALSE, sep = ",")

dat2 = read.csv("EC2_Expression_ProbesMatchedto_SF9495.csv")

all.equal(as.character(dat1$PROBE_ID), as.character(dat2$PROBE_ID))
# [1] TRUE

dim(dat2)
# [1] 47202    98

## Read in the new CorMatDNA function from Mike:

CorMatDNA = function(dat1, sampleindex1, simType, name1, WD1, dat2, sampleindex2, name2, WD2, subset, difference = c("dat1-dat2", "dat2-dat1", "reciprocal")) {
    timestamp()

    rtoz = function(x) {
        0.5 * log((1 + x) / (1 - x))
    }

    ztor = function(x) {
        (exp(2 * x) - 1) / (exp(2 * x) + 1)
    }

    if (class(dat1[, 1]) == "integer" | class(dat1[, 1]) == "numeric") {
        dat1 = dat1[order(as.numeric(as.character(dat1[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat1[, 1])))]
    } else {
        dat1 = dat1[order(as.character(dat1[, 1])), ]
        subset = subset[order(as.character(dat1[, 1]))]
    }

    if (class(dat2[, 1]) == "integer" | class(dat2[, 1]) == "numeric") {
        dat2 = dat2[order(as.numeric(as.character(dat2[, 1]))), ]
        subset = subset[order(as.numeric(as.character(dat2[, 1])))]
    } else {
        dat2 = dat2[order(as.character(dat2[, 1])), ]
        subset = subset[order(as.character(dat2[, 1]))]
    }

    if (all(dat1[, 1] == dat2[, 1]) == FALSE) {
        stop("dat1[,1] and dat2[,1] are not equivalent")
    }

    if (!is.null(subset)) {
        subset = as.logical(subset)
        datExpr1 = t(dat1[subset, sampleindex1])
        colnames(datExpr1) = as.character(dat1[subset, 1])
        datExpr2 = t(dat2[subset, sampleindex2])
        colnames(datExpr2) = as.character(dat2[subset, 1])
    } else {
        datExpr1 = t(dat1[, sampleindex1])
        colnames(datExpr1) = as.character(dat1[, 1])
        datExpr2 = t(dat2[, sampleindex2])
        colnames(datExpr2) = as.character(dat2[, 1])
    }

    if (simType == "Pearson") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "p", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "p", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Spearman") {
        print("Calculating corMat1...")
        corMat1 = cor(datExpr1, method = "s", use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = cor(datExpr2, method = "s", use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    if (simType == "Bicor") {
        print("Calculating corMat1...")
        corMat1 = bicor(datExpr1, use = "p")
        collectGarbage()
        diag(corMat1) = 0
        collectGarbage()
        print("Calculating corMat2...")
        corMat2 = bicor(datExpr2, use = "p")
        collectGarbage()
        diag(corMat2) = 0
        collectGarbage()
    }

    # print("Normalizing corMat1...")
    # corMat1=((ztor(scale(rtoz(corMat1))))+1)/2
    # corMat1=(ztor(rtoz(corMat1)/(1/sqrt(length(sampleindex1)-3)))+1)/2
    # diag(corMat1)=0
    # collectGarbage()

    if (min(corMat1, na.rm = T) < (-1) | max(corMat1, na.rm = T) > 1) {
        stop("Correlations in corMat1 outside of {-1,1}")
    }

    # print("Normalizing corMat2...")
    # corMat2=((ztor(scale(rtoz(corMat2))))+1)/2
    # corMat2=(ztor(rtoz(corMat2)/(1/sqrt(length(sampleindex2)-3)))+1)/2
    # diag(corMat2)=0
    # collectGarbage()

    if (min(corMat2, na.rm = T) < (-1) | max(corMat2, na.rm = T) > 1) {
        stop("Correlations in corMat2 outside of {-1,1}")
    }

    if (difference == "reciprocal") {
        simMat = corMat2 - corMat1
        filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        setwd(WD2)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
        simMat = corMat1 - corMat2
        filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        setwd(WD1)
        save(simMat, file = filenameout)
        rm(simMat, filenameout)
        collectGarbage()
    } else {
        if (difference == "dat2-dat1") {
            setwd(WD2)
            simMat = corMat2 - corMat1
            filenameout = paste("SD_", simType, "_", name2, "_minus_SD_", simType, "_", name1, sep = "")
        }

        if (difference == "dat1-dat2") {
            setwd(WD1)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")

            setwd(WD2)
            simMat = corMat1 - corMat2
            filenameout = paste("SD_", simType, "_", name1, "_minus_SD_", simType, "_", name2, sep = "")
        }

        collectGarbage()
        save(simMat, file = filenameout)
    }

    collectGarbage()

    timestamp()
} ## end of function


CorMatDNA(
    dat1 = dat1,
    sampleindex1 = c(7:75),
    name1 = "SF9495_PandG_30k",
    WD1 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    dat2 = dat2,
    sampleindex2 = c(7:98),
    name2 = "EC2_PandG_30k",
    WD2 = "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices",
    subset = !is.na(dat1$ProbeQuality) & (dat1$ProbeQuality == "P" | dat1$ProbeQuality == "G"),
    difference = "reciprocal",
    simType = "Bicor"
)


rm(list = ls())


## Now to build the consensus subtraction matrices:

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

load("SD_Bicor_SF9495_PandG_30k_minus_SD_Bicor_ACC1_PandG_30k")

simMat1 = simMat

rm(simMat)

dim(simMat1)
# [1] 30410 30410

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

load("SD_Bicor_SF9495_PandG_30k_minus_SD_Bicor_ACC2_PandG_30k")

simMat2 = simMat

rm(simMat)

dim(simMat2)
# [1] 30410 30410

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

load("SD_Bicor_SF9495_PandG_30k_minus_SD_Bicor_EC1_PandG_30k")

simMat3 = simMat

rm(simMat)

dim(simMat3)
# [1] 30410 30410

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices")

load("SD_Bicor_SF9495_PandG_30k_minus_SD_Bicor_EC2_PandG_30k")

simMat4 = simMat

rm(simMat)

dim(simMat4)
# [1] 30410 30410

## Just the 4 simMats in memory:

ls()
# [1] "simMat1" "simMat2" "simMat3" "simMat4"

## Since all of the subtraction matrices were constructed using the perfect and good probes from the array they are all the same size. I just need to check that the probe ID,s are the same:


all.equal(colnames(simMat1), colnames(simMat2))
# [1] TRUE
all.equal(colnames(simMat3), colnames(simMat4))
# [1] TRUE
all.equal(colnames(simMat1), colnames(simMat3))
# [1] TRUE
all.equal(colnames(simMat2), colnames(simMat4))
# [1] TRUE

all.equal(rownames(simMat1), rownames(simMat2))
# [1] TRUE
all.equal(rownames(simMat3), rownames(simMat4))
# [1] TRUE
all.equal(rownames(simMat1), rownames(simMat3))
# [1] TRUE
all.equal(rownames(simMat2), rownames(simMat4))
# [1] TRUE


## Generate a consensus of the 4 subtractions by taking the parallel minimum the 4 subtraction matrices:

conMat = pmin(simMat1, simMat2, simMat3, simMat4)

simMat = conMat

save(simMat, file = "Consensus_subtraction_SF9495_minus_ACC_and_EC")

## Remove the conMat and simMat from memory:
rm(conMat)
rm(simMat)



####################################### Building the consensus subtraction network ###########################################


## Move the consensus subtraction matrix to the network directory:
# file.copy("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices/Consensus_subtraction_SF9495_minus_ACC_and_EC", "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Consensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS_Modules/Consensus_subtraction_SF9495_minus_ACC_and_EC")
# [1] TRUE

## Using the latest version of FindModules.
# source("/fast-data/Shared/Code/KK/FindModules_0.90_KK.R")

# setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615")

# dat1=read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")

# dim(dat1)
# [1] 47202    75

# sum(!is.na(dat1$ProbeQuality)&(dat1$ProbeQuality=="P"|dat1$ProbeQuality=="G"))
# [1] 30410

# FindModules(
# projectname="Consensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS",
# expr=dat1,
# geneinfo=c(1:6),
# sampleindex=c(7:75),
# samplegroups=NULL,
# subset=!is.na(dat1$ProbeQuality)&(dat1$ProbeQuality=="P"|dat1$ProbeQuality=="G"),
# simMat="Consensus_subtraction_SF9495_minus_ACC_and_EC",
# saveSimMat=FALSE,
# simType="Bicor",
# beta=1,
# overlapType="None",
# TOtype="signed",
# TOdenom="min",
# MIestimator="mi.mm",
# MIdisc="equalfreq",
# signumType="rel",
# iterate=TRUE,
# signumvec=c(.9999,.999,.99,.98,.97,.96,.95,.94,.93,.92,.91,.90),
# minsizevec=c(5,6,7,8,9,10),
# signum=NULL,
# minSize=NULL,
# minMEcor=0.85,
# ZNCcut=2,
# calcSW=FALSE,
# loadTree=FALSE,
# writeKME=TRUE
# )

# rm(list=ls())




# .libPaths("/fast-data/R-3.1.1-packages")
# library(WGCNA)

# setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Consensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS_Modules")

# wd1=getwd()
## Get a list of directories from the wd:

# filelist1=list.dirs(wd1)

# for(i in 2:length(filelist1)){

#   setwd(filelist1[[i]])

#   wd1=getwd()

#   list1=list.files(wd1)

#   if(length(list1)==0){

# 		i=i+1

# 	} else {

#       kME1=grep("kME_table",list1)

#       datkme=read.csv(list1[[kME1]])

## A loop to prevent crashes when only one module is found:

#       if(length(grep("kME",colnames(datkme)))<4){

#           i=i+1

#       } else {

#           source("/fast-data/Shared/Code/SJS/GSEAfxsV1.r")

#           MyGSHGtest=MyGSHG(datkme1=datkme,kmecut1="topmodposbc",exclude="none",pvalcut1=NULL)
#           MyGSHGtest=MyGSHG(datkme1=datkme,kmecut1="topmodposfdr",exclude="none",pvalcut1=NULL)
#           MyGSHGtest=MyGSHG(datkme1=datkme,kmecut1="seed",exclude="none",pvalcut1=NULL)

#       }

#   }
# }




## Mike wanted me to rerun the GSEA on just the top 1000 probes in the turquoise and blue modules. I rewrote the GSEAFxr to allow you to assign a topx of genes. The modified version of the code is GSEAfxsV2.r and requires a topx to be input.


# source("/fast-data/Shared/Code/SJS/GSEAfxsV2.r")

# setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Consensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS_Modules/Bicor-None_signum0.621_minSize10_minMEcor0.85_30410")

# datkme=read.csv("kME_table_06-17-19.csv")

## Manually set top x to 1000 before running the function
# topx=1000

## Run the topx function for the top 1000 genes by BC and FDR.
# MyGSHGtest=MyGSHG(datkme1=datkme,kmecut1="topxBC",topx=1000,exclude="none",pvalcut1=NULL)
# MyGSHGtest=MyGSHG(datkme1=datkme,kmecut1="topxFDR",topx=1000,exclude="none",pvalcut1=NULL)

# rm(list=ls())

## Move the consensus subtraction matrix to the network directory:
file.copy("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Subtraction_Matrices/Consensus_subtraction_SF9495_minus_ACC_and_EC", "/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/XConsensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS_Modules/Consensus_subtraction_SF9495_minus_ACC_and_EC")
# [1] TRUE

## Using the latest version of FindModules. I created a version of FindModules in which relsignumvec=datx. I can then take the values from the previous network and build the exact same network:
source("/fast-data/Shared/Code/SJS/FindModules_0.90_KK_SetRelSigNum.r")

## Create a manual signum based on the previous network:
getSN = read.csv("/big-data/Sam/GBM_networks/Consensus_Subtraction_Sam_GBM1_Minus_Sam_Datasets_PandG_30K_SJS_Modules/Bicor-None_p1_Consensus_Subtraction_Sam_GBM1_Minus_Sam_Datasets_PandG_30K_SJS_30410_network_statistics_07-54-50.csv")

SN = as.list(getSN[3, 2:25])
name1 = as.list(getSN[2, 2:25])
name1 = unlist(name1)
name1 = as.character(name1)
name1 = as.numeric(name1)
name1 = name1 * 100
names(SN) = name1
SN = as.data.frame(SN)
setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615")
write.table(SN, file = "relsignumvec_OriginalNetwork.csv", col.names = TRUE, row.names = FALSE, sep = ",")

datx = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/relsignumvec_OriginalNetwork.csv")
colnames(datx) = gsub("X", "", colnames(datx))
colnames(datx) = gsub(".1", "", colnames(datx), fixed = TRUE)
colnames(datx) = paste(colnames(datx), "%", sep = "")

datx = as.list(datx)
datx = unlist(datx)

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615")

dat1 = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/SF9495_tumor/Renumbered_SF9495_ALL_69_Qnorm.csv")

dim(dat1)
# [1] 47202    75

sum(!is.na(dat1$ProbeQuality) & (dat1$ProbeQuality == "P" | dat1$ProbeQuality == "G"))
# [1] 30410

FindModules(
    projectname = "XConsensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS",
    expr = dat1,
    geneinfo = c(1:6),
    sampleindex = c(7:75),
    samplegroups = NULL,
    subset = !is.na(dat1$ProbeQuality) & (dat1$ProbeQuality == "P" | dat1$ProbeQuality == "G"),
    simMat = "Consensus_subtraction_SF9495_minus_ACC_and_EC",
    saveSimMat = FALSE,
    simType = "Bicor",
    beta = 1,
    overlapType = "None",
    TOtype = "signed",
    TOdenom = "min",
    MIestimator = "mi.mm",
    MIdisc = "equalfreq",
    signumType = "rel",
    iterate = TRUE,
    signumvec = c(.9999, .999, .99, .98, .97, .96, .95, .94, .93, .92, .91, .90),
    minsizevec = c(5, 6, 7, 8, 9, 10),
    signum = NULL,
    minSize = NULL,
    minMEcor = 0.85,
    ZNCcut = 2,
    calcSW = FALSE,
    loadTree = FALSE,
    writeKME = TRUE
)

rm(list = ls())

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/XConsensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS_Modules/Bicor-None_signum0.621_minSize10_minMEcor0.85_30410")

datkme = read.csv("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/XConsensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS_Modules/Bicor-None_signum0.621_minSize10_minMEcor0.85_30410/kME_table_07-09-38.csv")

.libPaths("/fast-data/R-3.1.1-packages-SJS")
library(WGCNA)

setwd("/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/XConsensus_Subtraction_SF9495_Minus_ACC_and_EC_PandG_30K_SJS_Modules")

wd1 = getwd()
## Get a list of directories from the wd:

filelist1 = list.dirs(wd1)

for (i in 2:length(filelist1)) {
    setwd(filelist1[[i]])

    wd1 = getwd()

    list1 = list.files(wd1)

    if (length(list1) == 0) {
        i = i + 1
    } else {
        kME1 = grep("kME_table", list1)

        datkme = read.csv(list1[[kME1]])

        ## A loop to prevent crashes when only one module is found:

        if (length(grep("kME", colnames(datkme))) < 4) {
            i = i + 1
        } else {
            source("/fast-data/Shared/Code/SJS/GSEAfxsV1.r")

            MyGSHGtest = MyGSHG(datkme1 = datkme, kmecut1 = "topmodposbc", exclude = "none", pvalcut1 = NULL)
            MyGSHGtest = MyGSHG(datkme1 = datkme, kmecut1 = "topmodposfdr", exclude = "none", pvalcut1 = NULL)
            MyGSHGtest = MyGSHG(datkme1 = datkme, kmecut1 = "seed", exclude = "none", pvalcut1 = NULL)
        }
    }
}
