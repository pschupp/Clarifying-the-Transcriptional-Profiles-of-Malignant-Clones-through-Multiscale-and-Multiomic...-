.libPaths("/home/shared/R/x86_64-pc-linux-gnu-library/4.0")
source("/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/downsample_fxns.R")
setwd("/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/SF10711/reads_per_sample/")
readsPerSample <- list.files()[grep("SF10711_9_", list.files())]
outDir <- "/home/rebecca/cancer_projects/astrocytoma/ampseq_downsample/SF10711/downsample_coverage/"
nResamples <- 1000
nThreads <- 15

############ SSH1 ############ 

projectName <- c("SSH1")
whichAmp <- c("SSH1_G1606T")
coverageIntervals <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 5000, 10000, 11000, 12000, 13000, 14000)

min(x$SSH1_G1606T)
# [1] 27571

downsampleSelfCov(projectName,
                  readsPerSample,
                  whichAmp,
                  coverageIntervals,
                  nResamples,
                  nThreads,
                  outDir)

############ TP53 ############ 

projectName <- c("TP53")
whichAmp <- c("TP53_G338T")
coverageIntervals <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 5000, 6000)

downsampleSelfCov(projectName,
                  readsPerSample,
                  whichAmp,
                  coverageIntervals,
                  nResamples,
                  nThreads,
                  outDir)

############ IDH1 ############ 

projectName <- c("IDH1")
whichAmp <- c("IDH1_G395A")
coverageIntervals <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 5000, 10000)

downsampleSelfCov(projectName,
                  readsPerSample,
                  whichAmp,
                  coverageIntervals,
                  nResamples,
                  nThreads,
                  outDir)
