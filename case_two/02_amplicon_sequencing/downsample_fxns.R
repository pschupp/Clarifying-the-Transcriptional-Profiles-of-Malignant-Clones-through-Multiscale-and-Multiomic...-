library(data.table)
library(dplyr)
library(WGCNA)
library(openblasctl)
library(doFuture)
library(doRNG)
library(flexiblas)
flexiblas_load_backend("OPENBLASPTHREAD")
flexiblas_switch(2)
options(future.globals.maxSize = +Inf)
openblas_set_num_threads(1)

downsampleCoverage <- function(projectName,
                               readsPerSample, ## Paths to reads per sample (one file per sample)
                               whichAmps,
                               coverageIntervals,
                               nResamples,
                               nThreads,
                               outDir) {
  registerDoFuture()
  plan(multicore, workers = nThreads)
  for (coverage in coverageIntervals) {
    cat("\n--- Starting coverage =", coverage,"---\n")
    downsampledCoverage <- list()
    for (sample in readsPerSample) {
      cat("\nStarting", sample, "\n")
      workingSample <- as.data.frame(fread(sample, nThread = nThreads))
      downsample <- foreach(seq(nResamples), .combine = rbind) %dorng% {
        ## Initialize downsample1 with first amp
        workingAmp <- workingSample[workingSample[, 1] %in% whichAmps[1],]
        replace <- FALSE
        if (nrow(workingAmp) < coverage)
          replace <- TRUE
        downsample1 <- sample_n(workingAmp, coverage, replace = replace) 
        ## Rest amps
        for (amp in whichAmps[2:length(whichAmps)]) {
          workingAmp <- workingSample[workingSample[, 1] %in% amp,]
          replace <- FALSE
          if(nrow(workingAmp) < coverage)
            replace <- TRUE
          downsample1 <- rbind(downsample1, sample_n(workingAmp, coverage, replace = replace)) ## Sample x reads per amplicon
        }
        ## Summarize number of alt and ref reads per amplicon:
        if(sum(downsample1$type %in% "alt") != 0) {
          downsample1 <- merge(as.data.frame(table(downsample1$gene[downsample1$type == "alt"])), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        } else {
          downsample1 <- merge(data.frame(whichAmps, NA), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        }
        downsample1[is.na(downsample1)] <- 0
        colnames(downsample1) <- c("Amplicon", "Alt", "Ref")
        return(downsample1)
      } # end downsamples
      downsampledCoverage[[sample]] <- downsample
      rm(downsample)
      collectGarbage()
    } # end samples
    assign(paste0("downsample_coverage_", coverage, "x"), downsampledCoverage)

  } ## // end coverage intervals
  save(list = ls()[grep("downsample_", ls())], file = paste0(outDir, "downsampleCoverage_", projectName, "_", nResamples, "_resamples.RData"))

} ## // end fxn

downsampleSelfCov <- function(projectName,
                              readsPerSample,
                              whichAmp,
                              coverageIntervals,
                              nResamples,
                              nThreads,
                              outDir) {
  registerDoFuture()
  plan(multicore, workers = nThreads)
  for (coverage in coverageIntervals) {
    cat("\n--- Starting coverage =", coverage,"---\n")
    downsampledCoverage <- list()
    for (sample in readsPerSample) {
      cat("\nStarting", sample, "\n")
      workingSample <- as.data.frame(fread(sample, nThread = nThreads))
      workingSample <- workingSample[workingSample[, 1] %in% whichAmp,]
      downsample <- foreach(seq(nResamples), .combine = cbind) %dorng% {
        idx <- sample(seq(nrow(workingSample)), coverage*2) 
        idx1 <- sample(idx, coverage)
        idx2 <- idx[!idx %in% idx1]
        downsample1 <- workingSample[idx1,]
        downsample2 <- workingSample[idx2,]
        if(sum(downsample1$type %in% "alt") != 0) {
          downsample1 <- merge(as.data.frame(table(downsample1$gene[downsample1$type == "alt"])), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        } else {
          downsample1 <- merge(data.frame(whichAmp, NA), as.data.frame(table(downsample1$gene[downsample1$type == "ref"])), by = 1, all = TRUE)
        }
        if(sum(downsample2$type %in% "alt") != 0) {
          downsample2 <- merge(as.data.frame(table(downsample2$gene[downsample2$type == "alt"])), as.data.frame(table(downsample2$gene[downsample2$type == "ref"])), by = 1, all = TRUE)
        } else {
          downsample2 <- merge(data.frame(whichAmp, NA), as.data.frame(table(downsample2$gene[downsample2$type == "ref"])), by = 1, all = TRUE)
        }
        downsample1[is.na(downsample1)] <- 0
        downsample2[is.na(downsample2)] <- 0
        colnames(downsample2) <- colnames(downsample1) <- c("Amplicon", "Alt", "Ref")
        return(rbind(downsample1, downsample2))
      } # end downsamples
      downsampledCoverage[[sample]] <- downsample
      rm(downsample)
      collectGarbage()
    } # end samples
    assign(paste0("downsample_self_coverage_", coverage, "x"), downsampledCoverage)
  } ## // end coverage intervals
  save(list = ls()[grep("downsample_self_coverage_", ls())], file = paste0(outDir, "downsampleSelfCov_", projectName, "_", nResamples, "_resamples.RData"))
} ## downsampleSelfCov()

downsampleSamples <- function(projectName,
                              readMat, ## First column = samples, 2 columns per amp corresponding to ref, alt
                              whichAmps,
                              nSamples,
                              nResamples,
                              outDir) {
  registerDoFuture()
  plan(multicore, workers = nThreads)
  for (n in nSamples) {
    cat("\n--- Starting samples =", n,"---\n")
    subsample <- foreach(seq(nResamples)) %dorng% {
      readMat1 <- readMat[sample(seq(nrow(readMat)), n),]
      freqs <- data.frame(Sample = readMat1[, 1])
      for (i in 1:length(whichAmps)) {
        workingAmp <- readMat1[, grep(whichAmps[i], colnames(readMat1))]
        freqs[, i + 1] <- workingAmp[, 2] / (workingAmp[, 2] + workingAmp[, 1])
      } ## // end amps
      colnames(freqs)[2:ncol(freqs)] <- whichAmps
      return(freqs)
    } ## // end resamples
    assign(paste0("downsample_", n, "_samples"), subsample)
  } ## // end n samples
  save(list = ls()[grep("downsample_", ls())], file = paste0(outDir, "downsampleSamples_", projectName, "_", nResamples, "_resamples.RData"))
} ## // end fxn

downsampleSampleskME <- function(projectName,
                                 readMat, ## First column = samples, 2 columns corresponding to ref, alt
                                 whichAmp,
                                 datExpr,
                                 truekME,
                                 nSamples,
                                 nResamples,
                                 outDir) {
      registerDoFuture()
      plan(multicore, workers = nThreads)
      kMECor <- list()
      for (n in nSamples) {
        cat("\n--- Starting samples =", n,"---\n")
        kME <- foreach(seq(nResamples)) %dorng% {
          readMat1 <- readMat[sample(seq(nrow(readMat)), n), c(1, grep(whichAmp, colnames(readMat)))]
          freq <- data.frame(Sample = readMat1[, 1],  readMat1[, 3] / (readMat1[, 3] + readMat1[, 2]))
          datExpr1 <- datExpr[, c(1, match(freq[, 1], colnames(datExpr)))]
          if(!identical(names(colnames(datExpr1)[-c(1)]), names(freq[, 1])))
            stop("\nSamples do not match\n")
          kME1 <- apply(datExpr1[, -c(1)], 1, function(x) cor(x, freq[, 2]))
          names(kME1) <- datExpr1[, 1]
          return(kME1)
        } ## // end resamples
        if(!identical(names(truekME), names(kME[[1]])))
          stop("\nGenes do not match\n")
        kMECor[[paste0(n, "_samples")]] <- mapply(function(x) cor(x, truekME, use = "complete.obs"), kME)
      } ## // end n samples
      save(kMECor, file = paste0(outDir, "downsampleSampleskME_correlation_", projectName, "_", nResamples, "_resamples.RData"))
    } ## // end fxn


downsampleCoverageAndSamples <- function(projectName,
                                         downsampledCoverage,
                                         whichAmps,
                                         nSamples,
                                         nResamples,
                                         outDir) {
  for (downsample in downsampledCoverage) {
    cat(paste("\n----- Starting", downsample, "-----\n"))
    workingDownsample <- eval(as.name(downsample))
    ## Subset each sample to amplicon 1 & 2
    workingAmp1 <- lapply(workingDownsample, function(x) x[x[, 1] %in% whichAmps[1],])
    workingAmp2 <- lapply(workingDownsample, function(x) x[x[, 1] %in% whichAmps[2],])
    ## Calculate VAF across all samples for amplicon 1 & 2 per resample
    freq1 <- mapply(function(resample) mapply(function(x) x$Alt[resample] / (x$Alt[resample] + x$Ref[resample]), workingAmp1), seq(nrow(workingAmp1[[1]])))
    freq2 <- mapply(function(resample) mapply(function(x) x$Alt[resample] / (x$Alt[resample] + x$Ref[resample]), workingAmp2), seq(nrow(workingAmp1[[1]])))
    freqCor <- vector(mode = "list", length(nSamples))
    names(freqCor) <- paste0(nSamples, "_samples")
    for (n in nSamples) {
      cat("\n--- Starting", n, "samples ---\n")
      for (i in 1:nResamples) {
        ## Fix coverage sampling, vary which samples
        freqCor1 <- c()
        for (j in 1:nResamples) {
          whichSamples <- sample(seq(nrow(freq1)), n)
          freqCor1 <- c(freqCor1, cor(freq1[whichSamples, i], freq2[whichSamples, i]))
        } ##  end j resamples of samples
        freqCor[[paste0(n, "_samples")]] <- c(freqCor[[paste0(n, "_samples")]], freqCor1)
      } ## // end i resamples of coverage
    } ## // end 
    assign(paste0("downsample_coverage_", gsub("downsample_", "", downsample)), freqCor)
  } ## // end downsamples
  save(list = ls(pattern = "downsample_coverage_"), file = paste0(outDir, "downsampleCoverageAndSamples_", projectName, "_", nResamples, "_resamples.RData"))
} ## // end fxn

downsampleCoverageAndSampleskME <- function(projectName,
                                            downsampledCoverage,
                                            whichAmp,
                                            datExpr,
                                            truekME,
                                            nSamples,
                                            nResamples,
                                            outDir) {
  registerDoFuture()
  plan(multicore, workers = nThreads)
  for (downsample in downsampledCoverage) {
    cat(paste("\n----- Starting", downsample, "-----\n"))
    workingDownsample <- eval(as.name(downsample))
    workingAmp <- lapply(workingDownsample, function(x) x[x[, 1] %in% whichAmp,])
    ## Calculate VAF across all samples per resample
    freq <- mapply(function(resample) mapply(function(x) x$Alt[resample] / (x$Alt[resample] + x$Ref[resample]), workingAmp), seq(nrow(workingAmp[[1]])))
    rownames(freq) <- gsub(".txt", "", gsub("sample_", "", rownames(freq)))
    kMECor <- list()
    ## Use all samples, vary coverage
    datExpr1 <- datExpr[, match(rownames(freq), colnames(datExpr))]
    freq1 <- freq[, 1:nResamples]
    if(!identical(colnames(datExpr1), rownames(freq)))
      stop("\nSamples do not match\n")
    kME <- apply(freq1, 2, function(x) apply(datExpr1, 1, function(y) cor(x, y)))
    if(!identical(datExpr[, 1], names(truekME)))
      stop("\nGenes do not match\n")
    kMECor1 <- apply(kME, 2, function(x) cor(x, truekME, use = "complete.obs"))
    kMECor[[paste0(nrow(freq), "_samples")]] <- kMECor1
    ## Vary samples and coverage
    for (n in nSamples) {
      cat("\n--- Starting", n, "samples ---\n")
      for (i in 1:nResamples) {
        ## Fix coverage sampling, vary which samples
        kMECor1 <- foreach(seq(nResamples), .combine = c) %dorng% {
          whichSamples <- rownames(freq)[sample(seq(nrow(freq)), n)]
          freq1 <- freq[match(whichSamples, rownames(freq)), i]
          datExpr1 <- datExpr[, match(whichSamples, colnames(datExpr))]
          if(!identical(colnames(datExpr1), names(freq1)))
            stop("\nSamples do not match\n")
          kME <- apply(datExpr1, 1, function(x) cor(x, freq1))
          if(!identical(datExpr[, 1], names(truekME)))
            stop("\nGenes do not match\n")
          return(cor(kME, truekME, use = "complete.obs"))
        } ##  end j resamples of samples
        kMECor[[paste0(n, "_samples")]] <- c(kMECor[[paste0(n, "_samples")]], kMECor1)
      } ## // end i resamples of coverage
    } ## // end 
    assign(paste0("downsample_coverage_kME_cor_", gsub("downsample_", "", downsample)), kMECor)
  } ## // end downsamples
  save(list = ls(pattern = "downsample_coverage_kME_cor_"), file = paste0(outDir, "downsampleCoverageAndSampleskME_correlation_", projectName, "_", nResamples, "_resamples.RData"))
} ## // end fxn
