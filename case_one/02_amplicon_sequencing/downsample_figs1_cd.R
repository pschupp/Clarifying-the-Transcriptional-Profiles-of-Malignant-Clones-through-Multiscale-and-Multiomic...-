.libPaths("/home/shared/R/x86_64-pc-linux-gnu-library/4.0")
source("/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/downsample_fxns.R")
setwd("/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/SF9495/reads_per_sample/")
readsPerSample <- list.files()
outDir <- "/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/SF9495/downsample_coverage/"
nResamples <- 1000
nThreads <- 12

############ TP53 ############ 

projectName <- c("TP53")
whichAmp <- c("TP53.L145P")
coverageIntervals <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

# min(x$TP53.L145P)
# [1] 2392

downsampleSelfCov(projectName,
                  readsPerSample,
                  whichAmp,
                  coverageIntervals,
                  nResamples,
                  nThreads,
                  outDir)

############ IDH1 ############ 

projectName <- c("IDH1")
whichAmp <- c("IDH1.R132H")
coverageIntervals <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000)

# min(x$IDH1.R132H)
# [1] 4440

downsampleSelfCov(projectName,
                  readsPerSample,
                  whichAmp,
                  coverageIntervals,
                  nResamples,
                  nThreads,
                  outDir)
