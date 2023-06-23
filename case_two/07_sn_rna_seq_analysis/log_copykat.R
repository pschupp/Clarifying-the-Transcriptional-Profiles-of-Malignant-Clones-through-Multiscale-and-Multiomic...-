library('copykat')
library('data.table')
library('ggplot2')
library('future')
library('future.apply')
library('qvalue')
library('dendextend')
library('RColorBrewer')
library('grid')
library('gridExtra')
library('limma')
library('ComplexHeatmap')
plan(multicore, workers=20)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
# get expression matrix and cluster assignments
# {{{
# read in expression matrix to determine wanted dimensions
expr = fread('~/@patrick/SF10711/sn.rna.seq/190809_SF013_oldham_july_1/08_expression_matrix/expr.ensg.counts.csv')
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

countThresh = future_apply(expr[,-c(1,2)], 1, function(x) length(which(x<0.5)) / (ncol(expr)-2))
threshInd = which(countThresh<0.9)
exprCD = expr[threshInd,-c(1,2)]
exprD = expr[threshInd,-c(1,2)]
exprSum = future_apply(exprCD, 1, sum)
exprGene = expr$Gene[threshInd]
exprENSG = expr$ENSG[threshInd]
rownames(exprD) = exprENSG
# }}}
# get cytobands
# {{{
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
# }}}
# get rna-seq clusters
# {{{
clusters = read.csv('~/@patrick/SF10711/cnv.analysis/sn.rna.seq/casper/clusters.csv')
nonmalignant = c(9, 11, 8, 4, 12, 6)
controlSampleIDs = clusters$X[clusters$x %in% nonmalignant]
# }}}
# run copykat
# {{{
source('/opt/copykat/R/copykat.R')
copyDat = copykat(rawmat = exprD, 
                    id.type = 'Ensemble',
                    cell.line = 'no',
                    ngene.chr = 20, 
                    LOW.DR = 0.05,
                    UP.DR = 0.10,
                    win.size = 25, 
                    norm.cell.names = controlSampleIDs,
                    KS.cut = 0.2, 
                    sam.name = "sf10711_sn_rna_w_ref", 
                    distance = "euclidean",
                    output.seg = FALSE,
                    plot.genes = T,
                    genome = 'hg20',
                    n.cores = 20)
# }}}
# summary by cytoband
copyDatBin = copyDat$CNAmat
copyDatSum = future_lapply(seq_len(nrow(cytoT)), function(i){
    temp = copyDatBin[(copyDatBin$chrom == cytoT$V1[i]) & (copyDatBin$chrompos > cytoT$V2[i]) & (copyDatBin$chrompos < cytoT$V3[i]), ]
    apply(temp[,-seq(1,3)], 2, mean)
})
copyDatSum = do.call(cbind, copyDatSum)
colnames(copyDatSum) = paste0(cytoT$V1, cytoT$V4)
write.csv(copyDatSum, file = 'copykat_cytoband_summary.csv')
