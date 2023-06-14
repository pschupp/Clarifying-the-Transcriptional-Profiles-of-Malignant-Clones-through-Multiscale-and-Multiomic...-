class=read.csv('/home/patrick/mount/big-data2/Patrick/july_seq/06.featureCounts.umi/cluster.assig.csv')
class[,2]=gsub('2', 'oligo', class[,2])
class[,2]=gsub('6', 'neuron', class[,2])
class[,2]=gsub('4', 'microglia', class[,2])
class[,2]=gsub('1', 'astrocyte', class[,2])
class[,2]=gsub('0', 'core.cancer1', class[,2])
class[,2]=gsub('3', 'prolif.outg', class[,2])
class[,2]=gsub('5', 'core.cancer2', class[,2])
write.table(class, row.names=F, col.names=F, file='infercnv.annotation.tsv', sep='\t', quote=F)

gtf.temp[,6]=make.names(gtf[,6], unique=T)
write.table(gtf.temp[,c(6,1,2,3)], row.names=F, col.names=F, sep='\t', file='infercnv.ordering.tsv', quote=F)

infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix="infercnv.proc.tsv", 
    annotations_file="infercnv.annotation.tsv", 
    delim="\t", 
    gene_order_file="infercnv.ordering.tsv", 
    ref_group_names=c('oligo', 'neuron', 'microglia', 'astrocyte'),
    out_dir=1"output_dir",  # dir is auto-created for storing outputs
    cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
    cluster_by_groups=T   # cluster
)
infercnv_obj = infercnv::run(infercnv_obj, cutoff=1,out_dir="output_dir", cluster_by_groups=F, denoise=T, HMM=F)
