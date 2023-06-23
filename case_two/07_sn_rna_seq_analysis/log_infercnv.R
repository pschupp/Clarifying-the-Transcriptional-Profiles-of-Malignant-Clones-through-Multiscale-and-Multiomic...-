library('infercnv')
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
write.table(exprD, row.names = F, sep = '\t', quote = F, file = 'infercnv.proc.tsv')
# }}}
# get rna-seq clusters
# {{{
clusters = read.csv('~/@patrick/SF10711/cnv.analysis/sn.rna.seq/casper/clusters.csv')
nonmalignant = c(9, 11, 8, 4, 12, 6)
controlSampleIDs = clusters$X[clusters$x %in% nonmalignant]
controlSample = data.frame(ID = colnames(expr), status = rep('malignant' length(colnames(expr))))
controlSample$status[controlSample$ID %in%controlSampleIDs] = 'normal'
write.table(controlSampleIDs, row.names = F, sep = '\t', quote = F, file = 'infercnv.annotation.tsv')
# }}}
# run infercnv 
write.table(gtf[,c(6,1,2,3)], row.names=F, col.names=F, sep='\t', file='infercnv.ordering.tsv', quote=F)

infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix="infercnv.proc.tsv", 
    annotations_file="infercnv.annotation.tsv", 
    delim="\t", 
    gene_order_file="infercnv.ordering.tsv", 
    ref_group_names=c('normal'),
    out_dir=1"output_dir",  # dir is auto-created for storing outputs
    cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
    cluster_by_groups=T   # cluster
)
infercnv_obj = infercnv::run(infercnv_obj, cutoff=1,out_dir="output_dir", cluster_by_groups=F, denoise=T, HMM=F)
