# Wheelplot figure generation
# dependencies and custom function
# {{{
library("circlize")
library("ComplexHeatmap")
library("stringr")
library("svglite")
library("extrafont")
library("gridBase")
library("gridtext")
vec2col = function(col1, col2, n, vec, col3 = NA) {
    if (class(vec) == "numeric") {
        naNumb = length(which(is.na(vec)))
        if (is.na(col3)) {
            mypal <- colorRampPalette(c(col2, col1))(n)
        } else {
            mypal <- colorRampPalette(c(col2, col3, col1))(n)
        }
        vec[which(vec == Inf)] = min(vec[-which(vec == Inf)])
        if (max(vec, na.rm = T) > 1.1) {
            names(mypal) = seq(0, 100, length = n)
        } else {
            names(mypal) = seq(0, 1, length = n)
        }
        out = rep("#5b5b5b", length(vec))
        i = 1
        for (ea in vec) {
            if (is.na(ea)) {
                i = i + 1
                next
            }
            ind = which.min(abs(as.numeric(names(mypal)) - ea))
            out[i] = mypal[ind]
            i = i + 1
        }
    }
    if (class(vec) == "factor") {
        out = rep(col1, length(vec))
        out[which(vec == "No" | vec == "F" | vec == "FALSE")] = col2
        out[which(vec == 0 | is.na(vec))] = "#5b5b5b"
    }
    out
}
# }}}
# read in data
# {{{
setwd("/home/patrick/@patrick/SF9495/wheelplot/inputs")
dat = read.csv("Complete list of unfiltered mutations.csv")
dat = read.csv("Shelton_Supplementary_Table_2_pgs_mod_new.csv")
expr = read.csv("Renumbered_SF9495_ALL_69_Qnorm.csv")
setwd("..")
expr = expr[, -c(1, 3, 4, 5, 6)]
expr = expr[-which(sapply(expr$Gene, nchar) == 0), ]
remove = rep(NA, nrow(expr))
i = j = 1
for (gene in expr$Gene[duplicated(expr$Gene)]) {
    if (nchar(gene) == 0) {
        next
    }
    ind1 = which(expr$Gene == gene)
    temp = expr[ind1, -1]
    ind2 = seq(1, length(ind1))[-which.max(apply(temp, 1, var))]
    j = length(ind2)
    remove[i:(i + j - 1)] = ind1[ind2]
    i = i + j
}
remove = remove[-which(is.na(remove))]
expr = expr[-remove, ]
rownames(expr) = expr$Gene
expr = expr[, -1]
exprMean = rowMeans(expr)
exprPerc = (order(exprMean, decreasing = T) / length(exprMean)) * 100
names(exprPerc) = names(exprMean)
dat = data.frame(dat, MeanExprPerc = exprPerc[match(dat$Gene, names(exprPerc))])
# }}}
# read in TCGA mutation data
# {{{
mut = read.csv("/mnt/clust/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/TCGA_LGG_data_030615/LGG_FINAL_ANALYSIS.aggregated.capture.tcga.uuid.curated.somatic.maf.csv")
bio1 = read.csv("/mnt/clust/fast-data/Shared/InputData/Sam_datasets_reprocessed_010615/TCGA_LGG_data_030615/Clinical/Biotab/nationwidechildrens.org_biospecimen_cqcf_lgg.txt", sep = "\t")
mut$TumorSample = unlist(lapply(strsplit(mut$Tumor_Sample_Barcode, "-"), function(x) paste(x[1:3], collapse = "-")))
bio1 = bio1[which(bio1$histological_type == "Astrocytoma"), ]
mut = mut[which(mut$TumorSample %in% bio1$bcr_patient_barcode), ]
mutOut = data.frame(hugoGene = mut$Hugo_Symbol, type = mut$Variant_Classification, chr = mut$Chromosome, start = mut$Start_position, end = mut$End_position, sample = mut$TumorSample, alt = mut$t_alt_count, ref = mut$t_ref_count)
mutOut = data.frame(mutOut, vaf = mutOut$alt / (mutOut$ref + mutOut$alt))
mutTot = aggregate(data.frame(mutOut$vaf), by = list(mutOut$hugoGene), FUN = "mean")
colnames(mutTot) = c("Gene", "vaf")
dat = data.frame(dat, MutFreq.TCGA.LGG.Frequency = mutTot$vaf[match(dat$Gene, mutTot$Gene)])
write.table(dat, col.names = TRUE, row.names = F, sep = ",", file = "SF9495.superfull.table.csv", quote = F)
# }}}
# prepare for use vep for hgvs naming and update hg38 coordinates
# {{{
for (ea in grep("Var.freq", colnames(dat))) {
    dat[which(is.na(dat[, ea])), ea] = 0
}
rownames(dat) = seq(1, nrow(dat))
dat = dat[-which(is.na(dat$Position)), ]
out = paste(dat2$Accession_unique, dat2$start - 1, gsub("-", "", dat2$REF), gsub("-", "", dat2$alt), sep = ":")
out = data.frame(chromosome = gsub("chr", "", dat$chr), start = dat$position, end = dat$position, allele = paste(dat$ref_allele, dat$alt_allele, sep = "/"), strand = rep("+", nrow(dat)))
out$end[grep("-", dat$ref_allele)] = out$end[grep("-", dat$ref_allele)] - 1
out$end = out$end + (nchar(dat$ref_allele) - 1)
write.table(out, col.names = F, row.names = F, sep = " ", file = "SF9495.positions.for.vep.txt", quote = F)
# bash
# /opt/vep_ensembl/ensembl-vep/vep	-i SF9495.positions.for.vep.txt  -o SF9495.vep.output.txt --cache \
# 									--dir_cache /home/shared/vep_cache/ --fasta /home/shared/vep_cache/fasta/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
# 									--hgvs --hgvsg --shift_hgvs 1 --symbol --refseq --fork 10 --pick --assembly GRCh38 --tab \
# 									--fields "Uploaded_variation,SYMBOL,Location,STRAND,REF_ALLELE,Allele,Gene,Feature,Consequence,SIFT,POLYPHEN,IMPACT,HGVSc,HGVSp,HGVSg"
vep = read.table("SF9495.vep.output.txt", row.names = NULL, header = T, sep = "\t")
vep$REF_ALLELE[grep(">", vep$HGVSc)] = gsub(">.*", "", gsub(".*[0-9]{1,}", "", vep$HGVSc[grep(">", vep$HGVSc)]))
vep$SYMBOL[8] = "CRIPAK"
vep$Feature[8] = "NM_175918.3"
vep$HGVSc[8] = "NM_175918.3:c.52_53insTGCCCATGTGGAGTGCCCGCCTGCTCACACA"
vep$HGVSp[8] = "CRIPAK:p.Cys18fsTer?"
vep = vep[match(gsub("\\..*", "", dat$Accession_unique), gsub("\\..*", "", vep$Uploaded_variation)), ]
dat = data.frame(vep$SYMBOL, gsub(":.*", "", vep$Location), gsub(".*:", "", vep$Location), vep$STRAND, vep$REF_ALLELE, vep$Allele, vep$Feature, vep$Consequence, vep$HGVSc, vep$HGVSp, dat$Var_freq_10, dat$Var_freq_34, dat$Var_freq_59, dat$Sanger_verified, dat$Amplicon_verified, dat$MeanExprPerc, dat$MutFreq.TCGA.LGG.Frequency)
colnames(dat) = c("Gene", "Chromosome", "Position", "Strand", "Reference_allele", "Alternate_allele", "Refseq_transcript_ID", "Consequence", "HGVSc", "HGVSp", "Section_10_VAF", "Section_34_VAF", "Section_59_VAF", "Verified_sanger", "Verified_ampseq", "Mean_expression_percentile", "TCGA_astrocytoma_vaf")
dat = dat[-which(apply(dat[, grep("Section.*VAF$", colnames(dat))], 1, function(x) length(which(x == 0))) == 3), ]
dat = dat[order(apply(dat[, grep("Section.*VAF$", colnames(dat))], 1, mean, na.rm = T), decreasing = T), ]
write.table(dat, file = "SF9495_wheelplot_table_full.tsv", sep = "\t", quote = F, row.names = F)
# }}}
# read in and tidy dataafter adjustment from vep
# {{{
setwd("~/@patrick/SF9495/wheelplot/")
dat = read.table("SF9495_wheelplot_table.tsv", sep = "\t", header = T, row.names = NULL)
dat = data.frame(dat, Label = paste(dat$Gene, gsub(".*:", "", dat$HGVSp), sep = ":"))
dat = dat[-grep("PHF8|MUC21|PCLO", dat$Gene), ]
dat$Label[which(dat$HGVSp == "-")] = paste(dat$Gene[which(dat$HGVSp == "-")], gsub(",.*", "", dat$Consequence[which(dat$HGVSp == "-")]), sep = ":")
dat$Label = gsub("_variant", "", dat$Label)
nMax = max(unlist(lapply(dat$Label, nchar)))
new = floor(ceiling((nrow(dat)) / .75) / 2)
oldRow = nrow(dat)
naDat = data.frame(matrix(NA, nrow = ceiling(ceiling(new) / 2), ncol = ncol(dat)))
colnames(naDat) = colnames(dat)
naDat$Label = make.names(naDat$Label, unique = T)
dat = data.frame(rbind(dat, naDat))
dat[, c(1, 2, 4, 5, 6, 7, 8, 9, 10, 14, 15, 18)] = lapply(dat[, c(1, 2, 4, 5, 6, 7, 8, 9, 10, 14, 15, 18)], as.factor)
dat$Label = factor(dat$Label, levels = dat$Label)
nameCol = rep("#5b5b5b", length(dat$Label))
nameCol[which(dat$Verified_ampseq == "Yes")] = "black"
nameCol[which(dat$Verified_ampseq == "No")] = "black"
names(nameCol) = dat$Label
names(nameCol) = gsub(":p\\.", " ", names(nameCol))
names(nameCol) = gsub(":", " ", names(nameCol))
names(nameCol) = gsub("_", " ", names(nameCol))
names(nameCol) = gsub("Leu", "L", names(nameCol))
names(nameCol) = gsub("Ile", "I", names(nameCol))
names(nameCol) = gsub("Pro", "P", names(nameCol))
names(nameCol) = gsub("Phe", "F", names(nameCol))
names(nameCol) = gsub("Arg", "R", names(nameCol))
names(nameCol) = gsub("His", "H", names(nameCol))
names(nameCol) = gsub("Ser", "S", names(nameCol))
names(nameCol) = gsub("Cys", "C", names(nameCol))
names(nameCol) = gsub("Val", "V", names(nameCol))
names(nameCol) = gsub("Thr", "T", names(nameCol))
names(nameCol) = gsub("Trp", "W", names(nameCol))
names(nameCol) = gsub("Ala", "A", names(nameCol))
names(nameCol) = gsub("Asp", "D", names(nameCol))
names(nameCol) = gsub("Asn", "N", names(nameCol))
names(nameCol) = gsub("Gly", "G", names(nameCol))
names(nameCol) = gsub("Gln", "N", names(nameCol))
names(nameCol) = gsub("Glu", "E", names(nameCol))
dat$Label = names(nameCol)
dat$Label = factor(dat$Label, levels = dat$Label)
# }}}
# create wheelplot
# {{{
setwd("~/@patrick/SF9495/figures")
svg("SF9495.wheelplot.svg", height = 10, width = 7.2)
plot.new()
circle_size = unit(1, "snpc")
pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
par(omi = gridOMI(), new = TRUE)
circos.clear()
circos.par(points.overflow.warning = F, track.height = 0.105, start.degree = 90, circle.margin = 1, track.margin = c(0.006, 0.006))
circos.initialize(as.factor(dat$Label), xlim = c(-1, 1))
circos.track(dat$Label, bg.col = vec2col("#4daf4a", "white", 100, dat$Section_10_VAF), ylim = c(-1, 1), bg.lwd = 1.5, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    colVec = nameCol[which(names(nameCol) == CELL_META$sector.index)]
    labelVec = CELL_META$sector.index
    circos.text(x = CELL_META$xcenter, y = ylim[1] + mm_y(6), col = colVec, labels = , labelVec, facing = "clockwise", niceFacing = TRUE, family = "NimbusSan", font = 2, cex = 1.2, adj = c(0, 0.5))
})
circos.track(dat$Label, bg.col = vec2col("#4daf4a", "white", 100, dat$Section_34_VAF), ylim = c(-1, 1), bg.lwd = 1.5)
circos.track(dat$Label, bg.col = vec2col("#4daf4a", "white", 100, dat$Section_59_VAF), ylim = c(-1, 1), bg.lwd = 1.5)
colTemp = vec2col("black", "white", 100, dat$Verified_sanger)
circos.track(dat$Label, bg.col = colTemp, ylim = c(-1, 1), bg.lwd = 1.5)
colTemp = vec2col("black", "white", 100, dat$Verified_ampseq)
circos.track(dat$Label, bg.col = colTemp, ylim = c(-1, 1), bg.lwd = 1.5)
circos.track(dat$Label, bg.col = vec2col("#377eb8", "white", 100, dat$TCGA_astrocytoma_vaf), ylim = c(-1, 1), force.ylim = T, bg.lwd = 1.5)
circos.track(dat$Label, bg.col = vec2col("#e31a1c", "white", 100, dat$Mean_expression_percentile), ylim = c(-1, 1), bg.lwd = 1.5)
rect(-2, -.025, -.015, 2, col = "white", border = NA)
textAdj = 0.13
textStart = 0.98
text(x = -0.70, y = textStart - 0 * textAdj, labels = str_pad("10", 32, "left"), family = "NimbusSan", cex = 1.1, col = "#4daf4a")
text(x = -0.86, y = textStart - 1 * textAdj, labels = str_pad("Exome VAF      34", 32, "left"), family = "NimbusSan", cex = 1.1, col = "#4daf4a")
text(x = -2.13, y = textStart - 1 * textAdj, labels = str_pad("{", 32, "left"), family = "NimbusSan", cex = 3.4, col = "#4daf4a")
text(x = -.71, y = textStart - 2 * textAdj, labels = str_pad("59", 32, "left"), family = "NimbusSan", cex = 1.1, col = "#4daf4a")
text(x = -.84, y = textStart - 3 * textAdj, labels = str_pad("Sanger verified", 32, "left"), family = "NimbusSan", cex = 1.1, col = "black")
text(x = -.87, y = textStart - 4 * textAdj, labels = str_pad("Amp-seq verified", 32, "left"), family = "NimbusSan", cex = 1.1, col = "black")
text(x = -1.03, y = textStart - 5 * textAdj, labels = str_pad("Gene-based VAF, TCGA astro.", 32, "left"), family = "NimbusSan", cex = 1.1, col = "#377eb8")
text(x = -1, y = textStart - 6 * textAdj, labels = str_pad("Mean expression percentile", 32, "left"), family = "NimbusSan", cex = 1.1, col = "#e31a1c")
legendTCGA = Legend(at = c(0, 1), col_fun = colorRamp2(c(1, 0), c("#377eb8", "white")), title = "Gene-based VAF\nTCGA astrocytomas", border = F, title_position = "topcenter", title_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"), labels_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"))
legendExome = Legend(at = c(0, 1), col_fun = colorRamp2(c(1, 0), c("#4daf4a", "white")), title = "Exome\nVAF", border = F, title_position = "topcenter", title_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"), labels_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"))
legendSanger = Legend(labels = c("Yes", "No", "No data"), title = "Sanger/Amp-seq\nverified", legend_gp = gpar(fill = c("black", "white", "#5b5b5b")), border = T, title_position = "topcenter", title_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"), labels_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"), row_gap = unit(2, "mm"))
legendExpr = Legend(at = c(0, 100), col_fun = colorRamp2(c(100, 0), c("#e31a1c", "white")), title = "Mean expr.\npercentile", border = F, title_position = "topcenter", title_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"), labels_gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "NimbusSan"))
upViewport()
draw(packLegend(legendExome, legendSanger, legendTCGA, legendExpr, direction = "horizontal", gap = unit(33, "mm")), y = unit(1, "npc") - circle_size, just = "top")
dev.off()
# }}}
