# 1) get copy number derived from FACETS
# {{{
source("/home/patrick/code/git/code_proc_SF10711/clone_pipeline/functions.R")
library("stringr")
library("gridExtra")
library("data.table")
library("RColorBrewer")
library("ggplot2")
cnv14 = read.csv("/mnt/bdata/@patrick/SF9495/exome_seq/facets/Glioma_14.cval.250.csv")
cnv39 = read.csv("/mnt/bdata/@patrick/SF9495/exome_seq/facets/Glioma_39.cval.250.csv")
cnv69 = read.csv("/mnt/bdata/@patrick/SF9495/exome_seq/facets/Glioma_69.cval.250.csv")
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
write.csv(cnvs, file = "exome_facets_cnvs.csv", quote = F, row.names = F)
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
# }}}
# 2) read in ChAMPS vectors of CNA
# {{{
# read in methylation data
cnvMeth = fread("~/@patrick/SF9495/methylation_array/LogR_ratio_tumor_vs_normal_brain.csv")
cnvMethOut = data.frame(matrix(nrow = 1, ncol = ncol(cnvMeth) + 1))
for (i in seq(1, nrow(cnvs_u))) {
    temp = cnvMeth[which(cnvMeth$chrom == cnvs_u$chrom[i] & cnvMeth$pos > cnvs_u$start[i] & cnvMeth$pos < cnvs_u$end[i]), ]
    temp2 = apply(temp[, -c(1, 2)], 2, median)
    cnvMethOut[i, ] = c(temp[1, c(1, 2)], temp[nrow(temp), 2], temp2)
}
colnames(cnvMethOut) = c("chrom", "start", "stop", colnames(cnvMeth)[-c(1, 2)])
rownames(cnvMethOut) = paste(cnvMethOut$chrom, cnvMethOut$start, sep = "_")
write.table(cnvMethOut, row.names = F, quote = F, sep = ",", file = "methylation_raw_from_ChAMPS_start_stop.csv")
cnvMethOut = cnvMethOut[, -c(1, 2, 3)]
relabel = read.csv("~/@patrick/SF9495/methylation_array/relabelled section 1to69.csv")
renameIndex = match(relabel$New, as.numeric(gsub(".*_", "", colnames(cnvMethOut))))
colnames(cnvMethOut) = relabel$Original[renameIndex[!is.na(renameIndex)]]
cnvMethOut = t(cnvMethOut)
# Ampseq is missing section 37
# cnvMethOut=cnvMethOut[-which(rownames(cnvMethOut)==37),]
# convert methylation data cellular frequencies of CN
cnvMethOutCor = data.frame(matrix(nrow = nrow(cnvMethOut), ncol = ncol(cnvMethOut)))
for (i in seq(1, ncol(cnvMethOut))) {
    methWork = cnvMethOut[which(rownames(cnvMethOut) %in% rownames(cnvColDf)), i]
    cnvMethOutCor[, i] = cnvMethOut[, i] * (mean(cnvColDf[cnvColDf[, i] != 0, i] / methWork[cnvColDf[, i] != 0]))
}
cnvMethOutCor[cnvMethOutCor < 0] = 0
rownames(cnvMethOutCor) = rownames(cnvMethOut)
colnames(cnvMethOutCor) = c("chr2p_del", "chr2q_del", "chr4p_del", "chr7_amp", "chr10p_amp", "chr12p_amp", "chr12q_del", "chr17p_del")
# write.table(cnvMethOutCor, row.names=T, quote=F, sep=',', file='methylation_champs_cnvs.csv')
# }}}
# 3) read in AmpSeq infered vectors of CNA
# not actually done, any suitable amplicons? No, no suitable amplicons except for TP53, otherwise ACCS on 11, GPR126, aka ADGRG6, on 6, PHF8 and  GPR173 on X
# {{{
# read in amplicon vaf table
amp = read.table("~/@patrick/SF9495/ampseq/Combined_variant_frequencies_run1_and_run2_with_readCounts.csv", sep = ",", header = T)
amp = amp[-(nrow(amp)), ]
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
setwd("~/@patrick/SF9495/clonality/pyclone/inputs")
amp = amp[which(amp$sample %in% rownames(cnvMethOut)), ]
names = paste("sf", str_pad(gsub(".*_", "", amp$sample), width = 3, pad = 0, side = "left"), sep = "")
# next need to take vector of TP53/2 - eigen(module1 from figure, manually select best genes)
amp = amp[, c(1, grep("freq", colnames(amp)))]
amp = amp[, c(1, grep("TEC|CTNNA|TP53|FAM193B|GDPD2|IAPP|IDH1|IAPP|IDH1", colnames(amp)))]
m1ME = apply(amp[, -c(1, grep("TP53", colnames(amp)))], 1, mean)
# }}}
# aligning all vectors into summary dataframe and writing it out as csv
# {{{
setwd("~/@patrick/SF9495/clonality/figures")
chr17ampVec = ((amp[, grep("TP53", colnames(amp))]) - m1ME) * 2
# create smooth version of dotplot
meth17 = smooth.spline(cnvMethOutCor$chr17p_del, df = 15)$y
out = data.frame(chr17ampVec, meth17)
colnames(out) = c("ampseq", "methylation")
# setwd('~/bdata/SF9495/figures')
# write.table(out, row.names=F, quote=F, sep=',', 'rna-seq_derived_cnvs_chr17.csv')
# 0.75
# }}}
# making correlation dotplots with copy number alteration as title and correlation R-value listed
# {{{
# create smoothed version of dotplot
cnvMethMin = cnvMethOutCor = apply(cnvMethOutCor, 2, function(x) smooth.spline(x, df = 15)$y)[grep("14|39|69", rownames(cnvMethOutCor)), ]
plot = data.frame(reshape2::melt(cnvColDf), reshape2::melt(cnvMethMin))[, -c(2, 4)]
oldham_theme = function() {
    theme(
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 40, color = "black", family = "NimbusSan", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "NimbusSan", margin = margin(r = 10)),
        axis.title.y = element_text(size = 40, face = "bold", family = "NimbusSan", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size = 40, face = "bold", family = "NimbusSan", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(t = -20, b = -40)),
        plot.subtitle = element_text(size = 40, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(t = 40, b = -40)),
        axis.line.x = element_line(size = 3),
        axis.line.y = element_line(size = 3),
        plot.margin = unit(c(0, 2, 1, 2), "lines"),
        legend.position = "right",
        legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size = 30, family = "NimbusSan")
    )
}
plot = plot[-grep("chr12", plot$Var2.1), ]
plot.bak = plot
plot = plot[-grep("chr17", plot$Var2.1), ]
WD = "~/@patrick/SF9495/figures/"
pdf(paste0(WD, "exome_meth.pdf"), width = 13, height = 13)
titleP = ggplot() +
    ggtitle("Concordance of exome and") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan"))
titleP2 = ggplot() +
    ggtitle("DNA methylation CNV calls") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan"))
titleP3 = ggplot() +
    ggtitle("N=3 sections, Pearson's r=0.92") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 40, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(b = -40)))
mainP = ggplot(plot, aes(x = value.1, y = value, color = Var2.1)) +
    geom_point(size = 10) +
    theme_classic() +
    oldham_theme() +
    scale_color_manual(labels = c("Chr2p Del", "Chr2q Del", "Chr4p Del", "Chr7 Gain", "Chr10p Gain", "Chr17p LOH"), values = brewer.pal(7, "Set1")[-6]) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.02, 1), breaks = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(1)) +
    labs(x = "Mean frequency (DNA Methylation)", y = "Mean frequency (Exome)") +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.spacing.x = unit(.5, "cm"),
        legend.key.size = unit(2, "cm"), legend.box = "vertical", legend.margin = margin()
    ) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE, ncol = 3))
grid.arrange(titleP, titleP2, titleP3, mainP, ncol = 1, heights = c(1.2 / 20, 1.2 / 20, 0.9 / 20, 16.7 / 20))
dev.off()
pdf(paste0(WD, "exome_meth_no_leg.pdf"), width = 13, height = 13)
titleP = ggplot() +
    ggtitle("Concordance of exome and") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan"))
titleP2 = ggplot() +
    ggtitle("DNA methylation CNV calls") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan"))
titleP3 = ggplot() +
    ggtitle("N=3 sections, Pearson's r=0.92") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 40, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(b = -40)))
mainP = ggplot(plot, aes(x = value.1, y = value, color = Var2.1)) +
    geom_point(size = 10) +
    theme_classic() +
    oldham_theme() +
    scale_color_manual(labels = c("Chr2p Del", "Chr2q Del", "Chr4p Del", "Chr7 Gain", "Chr10p Gain", "Chr17p LOH"), values = brewer.pal(7, "Set1")[-6]) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.02, 1), breaks = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(1)) +
    labs(x = "Mean frequency (DNA Methylation)", y = "Mean frequency (Exome)") +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        legend.spacing.x = unit(.5, "cm"),
        legend.key.size = unit(2, "cm"), legend.box = "vertical", legend.margin = margin()
    ) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE, ncol = 3))
grid.arrange(titleP, titleP2, titleP3, mainP, ncol = 1, heights = c(1.2 / 20, 1.2 / 20, 0.9 / 20, 16.7 / 20))
dev.off()

plot = plot.bak
out = data.frame(chr17ampVec[amp$sample %in% c(10, 34, 59)], plot$value.1[plot$Var2.1 == "chr17p_del"])
colnames(out) = c("ampseq", "exome")
# out$ampseq=out$ampseq/1.1
# cor(chr17ampVec, meth17)
pdf(paste0(WD, "chr17_amp_meth.pdf"), width = 13, height = 13)
titleP = ggplot() +
    ggtitle("Concordance of exome and") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan"))
titleP2 = ggplot() +
    ggtitle("amp-seq chr17p LOH calls") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(b = -40)))
titleP3 = ggplot() +
    ggtitle("N=3 sections, Pearson's r=0.99") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 40, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(b = -40)))
mainP = ggplot(out, aes(x = ampseq, y = exome)) +
    geom_point(size = 10) +
    theme_classic() +
    oldham_theme() +
    scale_color_manual(labels = c("Chr2p Del", "Chr2q Del", "Chr4p Del", "Chr7 Amp", "Chr10p Amp", "Chr17p LOH"), values = brewer.pal(7, "Set1")[-6]) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.29, 0.7), breaks = c(0.3, 0.7)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.29, 0.7), breaks = c(0.7)) +
    labs(y = "Mean frequency (Exome)", x = "Mean frequency (Amp-seq)")
grid.arrange(titleP, titleP2, titleP3, mainP, ncol = 1, heights = c(1.2 / 20, 1.2 / 20, 0.9 / 20, 16.7 / 20))
dev.off()

# }}}
# downsampling analysis
# {{{
WD = "~/@patrick/SF9495/figures/"
library("reshape2")
library("ggplot2")
library("gridExtra")
# TP53 RMSE plot
# {{{
base::load("/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF9495/downsample_coverage_freq_error/freqErrorCoverage_TP53_100_resamples.RData")
errPerCov = errPerCov[order(as.numeric(gsub("coverage_|x", "", names(errPerCov))))]
names(errPerCov) = gsub("coverage_", "", names(errPerCov))
df = melt(errPerCov)
colnames(df) <- c("RMSE", "Coverage")
df = df[!(df[, 2] %in% c("2000x", "3000x", "4000x")), ]
df$Coverage <- factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])
df$RMSE = as.numeric(df$RMSE)
tp53RMSE = ggplot(df, aes(x = Coverage, y = RMSE, group = Coverage)) +
    theme_classic() +
    oldham_theme() +
    geom_boxplot(fill = "white", color = "black", notch = TRUE, outlier.alpha = 0, lwd = 1.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.08), breaks = c(0, 0.08)) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(r = 4, l = 25, t = 0, b = 15),
        plot.title = element_text(margin = margin(t = 0, b = -40)),
        plot.subtitle = element_text(margin = margin(t = 60, b = 10))
    ) +
    labs(title = "Effect on TP53 VAF between\nfull and downsampled coverage", subtitle = "1000 resamples per coverage")
# }}}
# TP53 correlation plot
# {{{
base::load("/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF9495/downsample_coverage_freq_correlation/freqCorSelfCov_TP53_1000_resamples.RData", verbose = TRUE)
corPerCoverage = corPerCoverage[order(as.numeric(gsub("coverage_|x", "", names(corPerCoverage))))]
names(corPerCoverage) = paste0(gsub("coverage_", "", names(corPerCoverage)), "x")
df = melt(corPerCoverage)
colnames(df) = c("Correlation", "Coverage")
df$Coverage = factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])
df$Correlation = as.numeric(df$Correlation)
tp53Corr = ggplot(df, aes(x = Coverage, y = Correlation, group = Coverage)) +
    theme_classic() +
    oldham_theme() +
    geom_boxplot(fill = "white", color = "black", notch = TRUE, outlier.alpha = 0, lwd = 1.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(.7, 1), breaks = c(.7, .7, .7, .7, .7, .7, .7, .7, 1)) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "NimbusSan"),
        plot.margin = margin(l = 25, t = 40),
        axis.title.y = element_text(margin = margin(t = 0, r = 33, b = 0, l = 0))
    )
# }}}
pdf(paste0(WD, "tp53_downsample.pdf"), width = 13, height = 13)
grid.arrange(tp53RMSE, tp53Corr, ncol = 1)
dev.off()
# IDH1 RMSE plot
# {{{
base::load("/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF9495/downsample_coverage_freq_error/freqErrorCoverage_IDH1_100_resamples.RData")
errPerCov = errPerCov[order(as.numeric(gsub("coverage_|x", "", names(errPerCov))))]
names(errPerCov) = gsub("coverage_", "", names(errPerCov))
df = melt(errPerCov)
colnames(df) <- c("RMSE", "Coverage")
df = df[!(df[, 2] %in% c("3000x", "4000x")), ]
df$Coverage <- factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])
df$RMSE = as.numeric(df$RMSE)
idh1RMSE = ggplot(df, aes(x = Coverage, y = RMSE, group = Coverage)) +
    theme_classic() +
    oldham_theme() +
    geom_boxplot(fill = "white", color = "black", notch = TRUE, outlier.alpha = 0, lwd = 1.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.08), breaks = c(0, 0.08)) +
    theme(
        axis.title.x = element_blank(),
        plot.margin = margin(r = 4, l = 25, t = 0, b = 15),
        axis.text.x = element_blank(),
        plot.title = element_text(margin = margin(t = 0, b = -40)),
        plot.subtitle = element_text(margin = margin(t = 60, b = 10))
    ) +
    labs(title = "Effect on IDH1 VAF between\nfull and downsampled coverage", subtitle = "1000 resamples per coverage")
# }}}
# IDH1 correlation plot
# {{{
base::load("/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF9495/downsample_coverage_freq_correlation/freqCorSelfCov_IDH1_1000_resamples.RData", verbose = TRUE)
corPerCoverage = corPerCoverage[order(as.numeric(gsub("coverage_|x", "", names(corPerCoverage))))]
names(corPerCoverage) = paste0(gsub("coverage_", "", names(corPerCoverage)), "x")
df = melt(corPerCoverage)
colnames(df) = c("Correlation", "Coverage")
df$Coverage = factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])
df$Correlation = as.numeric(df$Correlation)
idh1Corr = ggplot(df, aes(x = Coverage, y = Correlation, group = Coverage)) +
    theme_classic() +
    oldham_theme() +
    geom_boxplot(fill = "white", color = "black", notch = TRUE, outlier.alpha = 0, lwd = 1.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.25, 1), breaks = c(0.25, 1)) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "NimbusSan"),
        plot.margin = margin(l = 25, t = 40),
        axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))
    )
# }}}
pdf(paste0(WD, "idh1_downsample.pdf"), width = 13, height = 13)
grid.arrange(idh1RMSE, idh1Corr, ncol = 1)
dev.off()



# IDH1 RMSE plot
# {{{
base::load("/mnt/bdata/rebecca/cancer_projects/UCSF_cases/ampseq_downsample/SF9495/downsample_coverage_freq_error/freqErrorCoverage_IDH1_100_resamples.RData")
errPerCov = errPerCov[order(as.numeric(gsub("coverage_|x", "", names(errPerCov))))]
names(errPerCov) = gsub("coverage_", "", names(errPerCov))
df = melt(errPerCov)
colnames(df) <- c("RMSE", "Coverage")
df$Coverage <- factor(df$Coverage, levels = unique(df$Coverage)[order(as.numeric(gsub("x", "", unique(df$Coverage))))])
df$RMSE = as.numeric(df$RMSE)
idh1RMSE = print(ggplot(df, aes(x = Coverage, y = RMSE, group = Coverage)) +
    theme_classic() +
    oldham_theme() +
    geom_boxplot(fill = "white", color = "black", notch = TRUE, outlier.alpha = 0, lwd = 1.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.08), breaks = c(0, 0.08)) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, family = "NimbusSan"),
        plot.title = element_text(margin = margin(t = 0, b = -40)),
        plot.subtitle = element_text(margin = margin(t = 60, b = 10)),
    ) +
    labs(title = "IDH1 VAF RMSE between\nfull and downsampled coverage", subtitle = "1000 resamples per coverage"))
# }}}
# TP53 correlation plot

# ddPCR plot
# {{{
library("reshape2")
library("ggplot2")
library("gridExtra")
WD = "~/@patrick/SF9495/figures/"
ddp1 = read.csv("~/cluster/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Digital_droplet_PCR/Final Analysis.csv")
ddp2 = read.csv("~/cluster/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Digital_droplet_PCR/Section Plate 2(Mut only).csv")
ddp = data.frame(rep1 = ddp1$FractionalAbundance[1:69], rep2 = ddp2$FractionalAbundance[1:69], rep1min = ddp1$PoissonFractionalAbundanceMin[1:69], rep2min = ddp2$PoissonFractionalAbundanceMin[1:69], rep1max = ddp1$PoissonFractionalAbundanceMax[1:69], rep2max = ddp2$PoissonFractionalAbundanceMax[1:69]) / 100
rownames(ddp) = gsub("Section ", "", ddp1$Sample)[1:69]
amp2 = read.table("/mnt/bdata/patrick/SF9495/ampseq/Combined_variant_frequencies_run1_and_run2_with_readCounts.csv", sep = ",", header = T)
amp2$IDH1.242.GtoA_UCI = amp2$IDH1.242.GtoA_freq + 2 * amp2$IDH1.242.GtoA_SE
amp2$IDH1.242.GtoA_LCI = amp2$IDH1.242.GtoA_freq - 2 * amp2$IDH1.242.GtoA_SE
amp2 = amp2[-(nrow(amp2)), ]
ddp = ddp[match(amp2$sample, rownames(ddp)), ]
ddpB = data.frame(sample = as.factor(rownames(ddp)), ddp[, c(1, 2)], amp2$IDH1.242.GtoA_freq, UCI = amp2$IDH1.242.GtoA_UCI, LCI = amp2$IDH1.242.GtoA_LCI)
colnames(ddpB)[1:4] = c("sample", "ddPCR rep. 1", "ddPCR rep. 2", "Amp-seq")
ddpM = reshape2::melt(ddpB[, c(1, 2, 3, 4)])
ddpM$sample = as.numeric(as.character(ddpM$sample))
ddpM$UCI = c(ddp$rep1max, ddp$rep2max, amp2$IDH1.242.GtoA_UCI)
ddpM$LCI = c(ddp$rep1min, ddp$rep2min, amp2$IDH1.242.GtoA_LCI)
pdf(paste0(WD, "ddpcr.pdf"), width = 16, height = 9.6)
titleP = ggplot() +
    ggtitle("IDH1 R132H") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan"))
mainP = ggplot(ddpM, aes(x = sample, y = value, group = variable, fill = variable, color = variable)) +
    theme_classic() +
    oldham_theme() +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, outline.type = "both", size = 0) +
    geom_line(size = 2) +
    labs(x = "Section ID", y = "VAF") +
    scale_x_continuous(expand = c(0, 0), limits = c(1, 81), breaks = c(1, 81), labels = c(1, 81)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.5), breaks = c(0.1, 0.5)) +
    scale_color_manual(values = c("#053061", "#2166ac", "#e41a1c")) +
    scale_fill_manual(values = c("#053061", "#2166ac", "#e41a1c")) +
    theme(
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.spacing.x = unit(1.0, "cm"),
        legend.key.size = unit(2, "cm"),
        plot.margin = margin(t = 0)
    )
grid.arrange(titleP, mainP, ncol = 1, heights = c(2 / 20, 18 / 20))
dev.off()
# }}}
# qPCR plot
# {{{
library("reshape2")
library("ggplot2")
library("gridExtra")
## The majority of the calculations for this analysis was done in excel. The csv file was then imported into R for me to make plots.
WD = "~/@patrick/SF9495/figures/"
CN = read.csv("~/cluster/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/Copy_number_analysis_QPCR/CopyNumberAnalysis.csv")
CN1 = data.frame(CN[, c(1, 8, 9)])
CN1 = data.frame(CN1[1:2], LCI = CN1[, 2] - 2 * CN1[, 3], UCI = CN1[, 2] + 2 * CN1[, 3])
CN1 = data.frame(CN1, Gene = "TP53")
colnames(CN1) = c("Section.No", "Copy.No", "LCI", "UCI", "Gene")
CN2 = data.frame(CN[, c(1, 16, 17)])
CN2 = data.frame(CN2[1:2], LCI = CN2[, 2] - 2 * CN2[, 3], UCI = CN2[, 2] + 2 * CN2[, 3])
CN2 = data.frame(CN2, Gene = "ACCS")
colnames(CN2) = c("Section.No", "Copy.No", "LCI", "UCI", "Gene")
CN = rbind(CN1, CN2)
CNSub1 = CN[c(1:69, 71:139), ]
CNBlood = CN[c(70, 140), ]
CNSub1$Section.No = as.character(CNSub1$Section.No)
CNBlood$Section.No = as.character(CNBlood$Section.No)
CNCombi = rbind(CNSub1, CNBlood)
CNSub1$Section.No = as.numeric(CNSub1$Section.No)
CNSub1$Copy.No = as.numeric(CNSub1$Copy.No)
pdf(paste0(WD, "qpcr_cnv_plot.pdf"), width = 16, height = 9.6)
titleP = ggplot() +
    ggtitle("No CNV identified by qPCR") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan"))
titleP2 = ggplot() +
    ggtitle("for ACCS and TP53") +
    geom_point() +
    theme_void() +
    theme(plot.title = element_text(size = 55, face = "bold", hjust = .5, family = "NimbusSan", margin = margin(b = 0)))
mainP = ggplot(CNSub1, aes(x = Section.No, y = Copy.No, color = Gene, fill = Gene, group = Gene)) +
    theme_classic() +
    oldham_theme() +
    scale_color_manual(values = c("#e41a1c", "#4daf4a")) +
    scale_fill_manual(values = c("#e41a1c", "#4daf4a")) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, outline.type = "both", size = 0) +
    geom_line(size = 2) +
    geom_hline(yintercept = 2, size = 2) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5), breaks = seq(0, 4)) +
    scale_x_continuous(expand = c(0, 0), limits = c(1, 81), breaks = c(1, 81), labels = c(1, 81)) +
    labs(y = "Copy number", x = "Section ID") +
    theme(
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.key.size = unit(2, "cm"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.spacing.x = unit(1.0, "cm"),
        plot.margin = margin(t = 0, l = 30)
    )
grid.arrange(titleP, titleP2, mainP, ncol = 1, heights = c(2.5 / 20, 2.5 / 20, 15 / 20))
dev.off()
# }}}
