# purpose of this is to redo casper
library('CaSpER')
library('GenomicRanges')
library('GenomeGraphs')
library('data.table')
library('biomaRt')
library('tidyr')
library('stringr')
library('openxlsx')
library('gdata')
library('preprocessCore')
library('data.table')
library('limma')
library('future')
library('future.apply')
library('circlize')
library('ComplexHeatmap')
source('/opt/CaSpER/R/plots.R')
source('/opt/Rphylip/Rphylip.R')
source('/opt/CaSpER/R/utility_functions.R')
plan(multicore, workers=20)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
# get expression matrix and cluster assignments
# {{{
# read in expression matrix to determine wanted dimensions
expr = data.frame(fread('~/@patrick/SF10711/sn.rna.seq/sanity/log_transcription_quotients.txt'))
expr1 =expr[,1]
expr = expr[,-1]
expr=data.table(expr1, apply(expr, 2, function(x) exp(x))*1E6)
# }}}
# harmonize gene names
# {{{
gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf', skip=5)
colnames(gtf)=c("chr", "db", "type", "start", "end", "score", "strand", "frame", "extra")
gtf=gtf[which(gtf$type=="gene"),]
gtf.ENSG=gsub(".*gene_id\\s\"\\\"*|\\..*", "", gtf$extra)
gtf.type=gsub(".*gene_type\\s\"\\\"*|\".*", "", gtf$extra)
gtf.name=gsub(".*gene_name\\s\"\\\"*|\".*", "", gtf$extra)
gtfOut=data.frame(gtf[,1], gtf.ENSG, gtf.type, gtf.name)
gtfOut[,4]= alias2SymbolTable(gtfOut[,4], species = 'Hs')

expr = data.frame(gtfOut[match(as.character(unlist(expr[,1])), gtfOut[,2]), 4], expr)
colnames(expr)[1:2] = c('Gene', 'ENSG')

exprC = fread('~/@patrick/SF10711/sn.rna.seq/190809_SF013_oldham_july_1/08_expression_matrix/expr.ensg.counts.csv')
exprC = exprC[match(expr$ENSG, exprC$Gene),]
exprC = data.frame(expr$Gene, exprC)
colnames(exprC)[1:2] = c('Gene', 'ENSG')

countThresh = future_apply(exprC[,-c(1,2)], 1, function(x) length(which(x<0.5)) / (ncol(exprC)-2))
threshInd = which(countThresh<0.9)
exprCD = exprC[threshInd,-c(1,2)]
exprD = expr[threshInd,-c(1,2)]
exprSum = future_apply(exprCD, 1, sum)
exprGene = expr$Gene[threshInd]
exprENSG = expr$ENSG[threshInd]
# }}}
# prep CaSPeR run
# {{{
thresh2 = is.na(exprGene)
exprD = exprD[!thresh2, ]
exprGene = exprGene[!thresh2]
exprENSG = exprENSG[!thresh2]
rownames(exprD) = exprENSG

centroL = read.table("~/@patrick/SF10711/cnv.analysis/sn.rna.seq/casper/cytoBand.txt", sep="\t")
centroTemp=centroL[,c(1,2,3,4)]
centroTemp[,1]=gsub("chr", "", centroTemp[,1])
centroTemp=centroTemp[which(nchar(centroTemp[,1])<=2),]
centroTemp[,4]= unlist(lapply(as.character(centroTemp[,4]), function(x) strsplit(x, "")[[1]][1]))
centroTemp=centroTemp[-grep("M", centroTemp[,1]),]
cytoT=data.frame(matrix(ncol=4, nrow=48))
i=1
for(chr in unique(centroTemp[,1])){
	for(cyto in c('p', 'q')){
		working=centroTemp[intersect(which(centroTemp[,1]==chr), which(centroTemp[,4]==cyto)),]
		cytoT[i,1]=chr
		cytoT[i,4]=cyto
		cytoT[i,2]=min(working[,2])+1
		cytoT[i,3]=max(working[,3])+1
		i=i+1
	}
}
cytoT=cytoT[order(as.numeric(cytoT[,1])),]
colnames(cytoT)=c("V1", "V2", "V3", "V4")

annoT =generateAnnotation(id_type="ensembl_gene_id", genes=exprENSG, ishg19=F, centroL, host="useast.ensembl.org")
lohL = readBAFExtractOutput(path='~/@patrick/SF10711/cnv.analysis/sn.rna.seq/casper', sequencing.type="single-cell", suffix="baf")
names(lohL) <- 'SF10711_p2'
lohNameMapping = data.frame(loh.name = rep("SF10711_p2", ncol(expr)-2), sample.name = colnames(expr)[-c(1,2)])
clusters = read.csv('~/@patrick/SF10711/cnv.analysis/sn.rna.seq/casper/clusters.csv')
nonmalignant = c(8, 7, 4, 10)
controlSampleIDs = clusters$X[clusters$x %in% nonmalignant]

exprD = exprD[match(annoT$Gene, rownames(exprD)),]
# exprGene = exprGene[match(annoT$Gene, exprGene)]

# }}}
# run CaSpER
# {{{
casperObj = CreateCasperObject(
                                raw.data = exprD,
                                annotation = annoT,
                                control.sample.ids = controlSampleIDs,
                                cytoband = cytoT,
                                loh.name.mapping = lohNameMapping,
                                cnv.scale = 3, 
                                loh.scale = 3,
                                method = 'iterative',
                                loh = lohL,
                                matrix.type = 'normalized',
                                sequencing.type = 'single-cell',
                                expr.cutoff = 0,
                                log.transformed = FALSE,
#                               centered.threshold = ,
                                window.length = 50,
                                length.iterations = 50,
                                vis.bound = 1,
                                genomeVersion = 'hg38'
                                )

casperRunObj = runCaSpER(casperObj, 
                            removeCentromere = T, 
                            cytoband = cytoT, 
                            method = 'iterative'
                        )
# }}}
# analyze CaSpER data 
setwd('~/@patrick/SF10711/cnv.analysis/sn.rna.seq/casper/analysis/nonmalig_reference')
# {{{
# Large-Scale CNV Summarization
finalChrMat = extractLargeScaleEvents(casperRunObj, thr=1) 
# }}}
# Segment based CNV Summarization
# {{{
gamma <- 6
all.segments <- do.call(rbind, lapply(casperRunObj, function(x) x@segments))
segment.summary <- extractSegmentSummary (casperRunObj)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh
loss.final <- loss[loss$count>gamma, ]
gain.final <- gain[gain$count>gamma, ]
loh.final <- loh[loh$count>gamma, ]
# }}}
# Gene based CNV Summarization
# {{{
all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), 
    IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(casperRunObj[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(geno.rna, ann.gr)
genes <- splitByOverlap(ann.gr, geno.rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(casperRunObj[[1]]@annotation.filt[,2])
all.samples <- unique(as.character(casperRunObj[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
# }}}
# get FM clone genesets and nonmalignant genesets
# {{{
clones = c('Clone.1', 'Clone.3', 'Clone.4', 'Clone.5')
fmGenes = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
clones=c('red', 'black', 'ivory', 'lightcyan')
fmGenes = lapply(clones, function(clone) fmGenes$Gene[which(fmGenes$'TopModPosFDR_0.12' == clone)])
names(fmGenes) = clones

fidGenes = lapply(list.files('/home/shared/genesets/fidelity_cell_types/hifi100', full.names=T), function(x) as.character(fread(x)[[1]]))
names(fidGenes) = gsub('_.*', '', list.files('/home/shared/genesets/fidelity_cell_types/hifi100'))
fidGenes= fidGenes[c(1, 10, 7, 9)]
# we can instead get the genes from the dataset
clones = c('Astrocyte', 'Oligodendrocyte', 'Microglia', 'Neuron')
bulkGenes = fread('~/@patrick/SF10711/integration_analysis/network_deconvolution/Bicor-None_signum0.438_minSize5_merge_ME_0.8_20246/kME_table_02-55-42.csv')
# clones=c('honeydew1', 'saddlebrown', 'coral1', 'salmon')
clones=c('honeydew1', 'tan', 'coral1', 'lavenderblush3')
bulkGenes = lapply(clones, function(clone) bulkGenes$Gene[which(bulkGenes$'TopModPosFDR_0.12' == clone)])
names(bulkGenes) = clones

set1 = c(fmGenes, bulkGenes)
# }}}

# get expression matrix and cluster assignments
# {{{
dist = as.numeric(unlist(fread('~/@patrick/SF10711/sn.rna.seq/sanity/cell_cell_distance_with_errorbar_s2n_gt_1.txt')))
# read in expression matrix to determine wanted dimensions
expr = data.frame(fread('~/@patrick/SF10711/sn.rna.seq/sanity/log_transcription_quotients.txt'))
expr1 =expr[,1]
expr = expr[,-1]
expr=data.table(expr1, apply(expr, 2, function(x) exp(x))*1E6)
# format into distance matrix
zeros = seq(1,ncol(expr)-1)
delta = seq((ncol(expr)-3), 1)
indeces = c(1)
for(i in seq(1, ncol(expr)-3)){
    if(i==1){
        indeces [i+1] = indeces[i] + delta[i]
    } else{
        indeces [i+1] = indeces[i] + delta[i] +1
    }
}
indeces[1]=0
indeces[length(indeces)+1]=indeces[length(indeces)]+1
distMat=lapply(seq_along(indeces[-length(indeces)]), function(i) {
    out=dist[(indeces[i]+1) : indeces[i+1]]
    if(zeros[i] != 809){
        out=unlist(c(rep(0, zeros[i]),out))
    } else{
        out=unlist(c(rep(0, zeros[i]-1), out[2]))
    }
    return(out)
})
distMat=data.frame(do.call(cbind, distMat))
distMat=cbind(distMat, rep(0, nrow(distMat)))
# now need to mirror distance matrix
for(i in seq_len(nrow(distMat))) {
    for(j in seq_len(ncol(distMat))) {
        distMat[i, j] = distMat[j, i]
    }
}
distDf = as.dist(distMat)
# }}}
# harmonize gene names
# {{{
gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf', skip=5)
colnames(gtf)=c("chr", "db", "type", "start", "end", "score", "strand", "frame", "extra")
gtf=gtf[which(gtf$type=="gene"),]
gtf.ENSG=gsub(".*gene_id\\s\"\\\"*|\\..*", "", gtf$extra)
gtf.type=gsub(".*gene_type\\s\"\\\"*|\".*", "", gtf$extra)
gtf.name=gsub(".*gene_name\\s\"\\\"*|\".*", "", gtf$extra)
gtfOut=data.frame(gtf[,1], gtf.ENSG, gtf.type, gtf.name)
gtfOut[,4]= alias2SymbolTable(gtfOut[,4], species = 'Hs')
fmGenes=lapply(fmGenes, alias2SymbolTable, species='Hs')

expr = data.frame(gtfOut[match(as.character(unlist(expr[,1])), gtfOut[,2]), 4], expr)
colnames(expr)[1:2] = c('Gene', 'ENSG')

exprC = fread('~/@patrick/SF10711/sn.rna.seq/190809_SF013_oldham_july_1/08_expression_matrix/expr.ensg.counts.csv')
exprC = exprC[match(expr$ENSG, exprC$Gene),]
exprC = data.frame(expr$Gene, exprC)
colnames(exprC)[1:2] = c('Gene', 'ENSG')

countThresh = future_apply(exprC[,-c(1,2)], 1, function(x) length(which(x<0.5)) / (ncol(exprC)-2))
threshInd = which(countThresh<0.9)
exprCD = exprC[threshInd,-c(1,2)]
exprD = expr[threshInd,-c(1,2)]
exprSum = future_apply(exprCD, 1, sum)
exprGene = expr$Gene[threshInd]
# }}}
# get cluster definitions from hierarchical clustering
# {{{
attr(distDf, 'Labels') =paste0('x', seq(1, length(attr(distDf, 'Labels'))))
distHclust = hclust(distDf, method = 'ward.D')
dendP = as.dendrogram(distHclust)
clust = cutree(distHclust, k=10)
names(clust) = colnames(expr)[-c(1,2)]
# }}}
# get single-nucleus-amplicon information
# {{{
snAmp=read.csv('~/@patrick/SF10711/sn.amp.seq/old/21_02_25_malignancy.adjusted.stringent.csv')
snAmp$id = gsub('(N[0-9]*)\\.(.*)', '\\2\\.\\1', snAmp$id)
snAmp=snAmp[match(colnames(exprD), snAmp$id),]
snAmp[(clust==10 & snAmp[,2]==1),2]='0'
snAmp[(clust==8 & snAmp[,2]==1),2]='0'
snAmp[(clust==4 & snAmp[,2]==1),2]='0'
snAmp[(clust==7 & snAmp[,2]==1),2]='0'
snAmp[(clust==3 & snAmp[,2]==0),2]='1'
snAmp[(clust==5 & snAmp[,2]==0),2]='1'
snAmp[(clust==1 & snAmp[,2]==0),2]='1'
snAmp[(clust==6 & snAmp[,2]==0),2]='1'
snAmp[(clust==9 & snAmp[,2]==0),2]='1'
snAmp[(clust==2 & snAmp[,2]==0),2]='1'
clust2 = cutree(distHclust, k=12)
snAmp[(clust2==6 & snAmp[,2]!='UNK'),2]='0'
snAmp$malignant = gsub('0', 'Nonmalignant', snAmp$malignant)
snAmp$malignant = gsub('1', 'Malignant', snAmp$malignant)
snAmp$malignant = gsub('UNK', 'No data', snAmp$malignant)
# }}}


# plots
# {{{
# heatmap
obj <- casperRunObj[[9]]
plotHeatmap10x(object=obj, fileName="heatmap.png",cnv.scale= 1, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = F)
plotLargeScaleEvent (object=obj, fileName="large.scale.events.png")
results <- extractMUAndCooccurence (finalChrMat, lohL, lohNameMapping)
resultsFil = list(results[[1]][which(results[[1]]$Pval<1E-13), ])
## visualize mutual exclusive and co-occurent events
plotMUAndCooccurence (resultsFil)
plotSCellCNVTree(finalChrMat, 'SF10711_p2', fileName = 'treeTry')

# visualize raw data
# fineCNV = casperRunObj[[9]]@control.normalized.noiseRemoved[[1]]
# # aggregate by cytoband
# fineCNVcyto = future_apply(fineCNV, 2, function(x) aggregate(x, by = list(annoT$band), FUN=mean)[,2])
# cytoAgg = aggregate(x, by = list(annoT$band), FUN=mean)
# rownames(fineCNVcyto) = cytoAgg[,1]
# fineCNVcyto = cbind(fineCNVcyto)
# fineCNVcyto = fineCNVcyto[match(unique(annoT$band), rownames(fineCNVcyto)), ]
# fineCNVcyto=(fineCNVcyto-1)*100+2
# fineCNVcyto[fineCNVcyto>1 & fineCNVcyto<2] = 2
# fineCNVcyto[fineCNVcyto>2 & fineCNVcyto<3] = 2
# pdf('try.pdf', width=15)
# Heatmap(t(fineCNVcyto),
#         show_column_names = FALSE,
#         show_row_names = FALSE,
#         col = colorRamp2(breaks = c(0, 2, 4),
#                             colors = c('#2166ac', '#f7f7f7', '#b2182b')
#                         ),
#         top_annotation = ha, 
#         left_annotation = rha, 
#         column_split = chroms,
#         row_split = 10,
#         cluster_columns = FALSE,
#         clustering_distance_rows = 'spearman',
#         clustering_method_rows = 'ward.D')
# dev.off()
# }}}

samps <- obj@large.scale.cnv.events
chrs <- as.vector(sapply(1:22, function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
chrMat <- matrix(0, ncol = 44, nrow = length(rownames(samps)))
colnames(chrMat) <- chrs
rownames(chrMat) <- rownames(samps)


for (x in 1:dim(samps)[1]) {
    
    chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleAmp[x], collapse = " ")), split = " ")))
    chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- 1
    
    chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleDel[x], collapse = " ")), split = " ")))
    chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- (-1)
    
}

plot.data <- melt(chrMat)
plot.data$value2 <- "neutral"
plot.data$value2[plot.data$value > 0] <- "amplification"
plot.data$value2[plot.data$value < 0] <- "deletion"
plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
plotW = reshape(plot.data[,c(1,2,4)], idvar = 'Var1', timevar ='Var2', direction = 'wide')
plotW2 = reshape(plot.data[,c(1,2,3)], idvar = 'Var1', timevar ='Var2', direction = 'wide')
rownames(plotW) = plotW[,1]
plotW = plotW[, -1]
colnames(plotW) = gsub('value2\\.', '', colnames(plotW))
rownames(plotW2) = plotW2[,1]
plotW2 = plotW2[, -1]
plotW2D = hclust(dist(plotW2, method = 'euclidean'), method = 'ward.D2')
cutreeCasp = cutree(plotW2D, k=2)
write.csv(cutreeCasp, file = 'CaSpER_clusters.csv')
colnames(plotW2) = gsub('value2\\.', '', colnames(plotW2))
plotCols = c('amplification' = '#b2182b', 'neutral' = '#f7f7f7', 'deletion' = '#2166ac')

chroms = annoT$Chr[match(rownames(fineCNVcyto), annoT$band)]
top_annotation = HeatmapAnnotation(chrom = anno_block(gp = gpar(fill = 2:4),
    labels = chroms, 
    labels_gp = gpar(col = "white", fontsize = 10)))
chroms = factor(chroms, levels = seq(1,21))
colsP = rep(c('grey','lightgrey'),length=19)
names(colsP) = unique(chroms)
ha = HeatmapAnnotation(chrom = chroms,
    col = list(chrom = colsP ),
    show_annotation_name = TRUE,
    show_legend = F)
clustCols = as.factor(c('#701d02', '#1b9e77', '#e7298a', '#4daf4a', '#4daf4a', '#999999', '#dd95e8', '#e6ab02', '#a65729', '#ff7f00'))
rha = rowAnnotation(
    'Amplicon-seq' = snAmp$malignant,
    'snRNA-seq clusters' = as.factor(clust),
    col = list('Amplicon-seq' = c('Nonmalignant' = '#2166ac', 'Malignant' = '#b2182b', 'No data' = '#999999'),
                'snRNA-seq clusters' = clustCols)
)
pdf('try.pdf')
Heatmap(plotW,
        col = plotCols,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        cluster_rows = plotW2D,
        left_annotation = rha)
#         clustering_distance_rows = 'binary',
#         clustering_method_rows = 'ward.D')
dev.off()
