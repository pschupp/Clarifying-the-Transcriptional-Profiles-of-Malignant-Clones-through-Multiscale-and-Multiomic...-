# dependencies
# {{{
library('data.table')
library('future')
library('future.apply')
library('data.table')
library('ggplot2')
library('ggnetwork')
library('STRINGdb')
library('png')
library('igraph')
library('RCurl')
library('ggnetwork')
library('network')
library('intergraph')
library('igraph')
library('gridExtra')
library('flashClust')
library('Hmisc')
library('RColorBrewer')
source('~/code/git/GSEA_generic/enrichment_functions.R')
source('/opt/stringdb/rstring.R')
source('/opt/stringdb/r_utils.R')
source("GSEAfxsV3.r")
plan(multicore, workers=20)
options(future.globals.maxSize = as.numeric('+Inf'))
options(future.rng.onMisuse='ignore')
broadSets=getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
# }}}
# method functions from package STRING DB rewritten
# {{{
map = function(my_data_frame,my_data_frame_id_col_names,takeFirst=TRUE, removeUnmappedRows=FALSE, quiet=FALSE){
        aliasDf2=get_aliases(takeFirst)
        tempDf = multi_map_df(my_data_frame, aliasDf2, my_data_frame_id_col_names, "alias", "STRING_id")
        naDf = subset(tempDf, is.na(STRING_id))
        if(nrow(naDf) > 0 & !quiet) cat(paste("Warning:  we couldn't map to STRING ", as.integer((nrow(naDf)/nrow(tempDf))*100), "% of your identifiers" , sep=""))
        if(removeUnmappedRows) tempDf = subset(tempDf, !is.na(STRING_id))
        return(tempDf)
      }
get_aliases = function(takeFirst=TRUE){
    temp = downloadAbsentFile(paste(protocol, "://stringdb-static.org/download/protein.aliases.v", file_version, "/", species, ".protein.aliases.v", file_version, ".txt.gz", sep="")) #, oD=input_directory)        
    if(!takeFirst){ 
      aliases_type <<- "all"
    } else {
      aliases_type <<- "take_first"
    }
    proteins <<- get_proteins()
    aliasDf <- read.table(temp, skip=1, sep = "\t", header=FALSE, quote="", stringsAsFactors=FALSE, fill = TRUE)
    colnames(aliasDf) <- c("STRING_id", "alias", "sources")
    aliasDf = subset(aliasDf, select=c("STRING_id", "alias"))
    pr1=data.frame(STRING_id=proteins$protein_external_id, alias=proteins$preferred_name, stringsAsFactors=FALSE)        
    pr2=data.frame(STRING_id=proteins$protein_external_id, alias=proteins$protein_external_id, stringsAsFactors=FALSE)
    pr3=data.frame(STRING_id=proteins$protein_external_id, alias=unlist(strsplit(proteins$protein_external_id, "\\."))[seq(from=2, to=2*nrow(proteins), by=2)], stringsAsFactors=FALSE)
    if(takeFirst){aliasDf = subset(aliasDf, !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(proteins$preferred_name)) & 
                                     !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(proteins$protein_external_id))  &
                                      !(toupper(iconv(alias, "WINDOWS-1252", "UTF-8")) %in% toupper(unlist(strsplit(proteins$protein_external_id, "\\."))[seq(from=2, to=2*nrow(proteins), by=2)])) )
    }
    aliasDf2=rbind(pr1,pr2,pr3, aliasDf)
    aliases_tf <<- aliasDf2
    return(aliasDf2)
}
get_proteins = function(){
if(nrow(proteins)==0){
    temp = downloadAbsentFile(paste(protocol, "://stringdb-static.org/download/protein.info.v", file_version, "/", species, ".protein.info.v", file_version, ".txt.gz", sep=""))# , oD=input_directory)
    if (version %in% c("11.0", "11.0b")) {
        proteinsDf <- read.table(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE, quote="")
        proteinsDf2 = subset(proteinsDf, select=c("protein_external_id",  "preferred_name", "protein_size", "annotation"))
    } else {
       proteinsDf <- read.table(temp, sep = "\t", skip=1, header=FALSE, stringsAsFactors=FALSE, fill = TRUE, quote="")
       colnames(proteinsDf) <- c("protein_external_id",  "preferred_name", "protein_size", "annotation")
       proteinsDf2 = subset(proteinsDf, select=c("protein_external_id",  "preferred_name", "protein_size", "annotation"))
    }
    proteins <<- proteinsDf2
  }
  return(proteins)
}
# }}}
# method functions from package STRING DB rewritten
# {{{
plot_network = function(string_ids, payload_id=NULL, required_score=NULL, add_link=FALSE, network_flavor='physical', add_summary=TRUE) {
    if (version %in% c("11.0")) {
        print("Parameter add_link not available in version 11.0 (please use 11.0b or later)")
        add_link = FALSE
    }
    if(is.null(required_score) ) required_score = score_threshold
    img = get_png(string_ids, payload_id=payload_id, required_score=required_score, network_flavor=network_flavor)
    if(!is.null(img)){
        plot(1:(dim(img)[2]), type='n', xaxt='n', yaxt='n', xlab="", ylab="", ylim=c(1,dim(img)[1]), xlim=c(1,(dim(img)[2])), asp = 1 )
        if(add_summary) mtext(get_summary(string_ids, required_score), cex = 0.7)
        if(add_link) mtext(get_link(string_ids, payload_id=payload_id, required_score=required_score), cex = 0.7, side=1)
        rasterImage(img, 1, 1, (dim(img)[2]), dim(img)[1])
    } 
}
get_png = function(string_ids, required_score=NULL, network_flavor="physical", file=NULL, payload_id=NULL){
    if(length(string_ids) > 2000) {
        cat("ERROR: We do not support lists with more than 2000 genes.\nPlease reduce the size of your input and rerun the analysis. \t")
        stop()
    }
    if(is.null(required_score) ) required_score = score_threshold
    string_ids = unique(string_ids)
    string_ids = string_ids[!is.na(string_ids)]
    urlStr = paste(stable_url, "/api/image/network", sep="")
    identifiers=""
    for(id in string_ids ){ identifiers = paste(identifiers, id, sep="%0d")}
    params = list(required_score=900, required_score=required_score, network_flavor='physical', identifiers=identifiers, species=species, caller_identity='STRINGdb-package')
    if(!is.null(payload_id)) params["internal_payload_id"]= payload_id
    img <- readPNG(postFormSmart(  urlStr, .params=params) )
    if(!is.null(file))  writePNG(img,  file)
    return(img)
}
get_summary = function(string_ids, required_score){
    string_ids = unique(string_ids)
    string_ids = string_ids[!is.na(string_ids)]
    enrichment = ppi_enrichment(string_ids, required_score)
    summaryText = paste("proteins: ", length(string_ids), '\n', 
                            "interactions: ", enrichment$edges, '\n',
                            "expected interactions: ", enrichment$lambda,' (',
                            "p-value: ",enrichment$enrichment, ')\n', sep=""
    )
    return(summaryText)
}
ppi_enrichment = function(string_ids, required_score=NULL){
    if (is.null(required_score)) required_score = score_threshold
    urlStr <- paste(stable_url, "/api/tsv-no-header/ppi_enrichment", sep="")
    string_ids = unique(string_ids)
    string_ids = string_ids[!is.na(string_ids)]
    identifiers = ""
    for(id in string_ids) {
      identifiers = paste(identifiers, id, sep="%0d")
    }
    if(length(backgroundV)==0) {
       params = list(species=species, identifiers=identifiers, required_score=required_score)
    }else {
        background = ""
        for(id in backgroundV) {
            background = paste(background, id, "%0d", sep="")
        }
        params = list(species=species, identifiers=identifiers, background_string_identifiers=background, required_score=required_score)
    } 
    tempDfv=postFormSmart(urlStr, .params=params)
    answer = read.table(text=tempDfv, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE, header=FALSE)
    result = list(enrichment = answer$V6, edges=answer$V2, lambda=answer$V5)
    return(result)
}
# }}}
# fisher transform functions
# {{{
rtoz = function(x) {
    0.5 * log((1 + x) / (1 - x))
}
calc.SD.z= function(no.samples, no.datasets) {
    1 / sqrt(no.samples - 3 * no.datasets)
}
ztor = function(x) {
    (exp(2 * x) - 1) / (exp(2 * x) + 1)
}
fisher.transform = function(x) {
    Mean.z = mean(rtoz(x))
    SD.z = calc.SD.z(no.samples = 65 + 85, no.datasets = 2)
    z.Fid = Mean.z / SD.z
    P.value = 2 * pnorm(-abs(z.Fid))
    Mean.r = ztor(Mean.z)
    P.value
}
fisher.transformtocor = function(x) {
    Mean.z = mean(rtoz(x))
    SD.z = calc.SD.z(no.samples = 65 + 85, no.datasets = 2)
    z.Fid = Mean.z / SD.z
    P.value = 2 * pnorm(-abs(z.Fid))
    Mean.r = ztor(Mean.z)
    Mean.r
}
# }}}
#####################################################################################################################################################
cor_method = 'correlation'
method = 'fisher'
pValAll=c()
# calculate correlation
# normalized SF9495 expression matrix - rnaScale9495
# {{{
rna=data.frame(fread('~/@patrick/SF9495/rna_array/Renumbered_SF9495_ALL_69_Qnorm.csv'))
meanExpr=future_apply(rna[,-seq(1,6)],1, mean)
rna=rna[order(meanExpr, decreasing=T),]
rna=rna[!(duplicated(rna$Gene)),]
rna=rna[-grep('^MIR', rna$Gene),]
rna=rna[-grep('^NCRNA', rna$Gene),]
rna=rna[-grep('^LOC', rna$Gene),]
rna=rna[-grep('^SNOR', rna$Gene),]
rna=rna[-grep('^DKFZ', rna$Gene),]
rna=rna[-grep('^KIAA', rna$Gene),]
rna=rna[-grep('^FLJ', rna$Gene),]
rnaScale9495=rna[,-c(seq(1,6),69,70)]
colnames(rnaScale9495)=gsub('SF9495_', '', colnames(rnaScale9495))
barOut = data.frame(fread('~/@patrick/SF9495/figures/sf9495_cellular_abundance.csv'))
rnaScale9495 = rnaScale9495[, c(1, which(as.numeric(colnames(rnaScale9495)) %in% barOut[,1])+1)]
rnaScale9495=data.frame(Gene=rna$Gene, rnaScale9495)
# }}}
# clone 1 abundance SF9495 - barOutIn9495
# {{{
barOutIn=data.frame(barOut[,-c(1,2,3,4,5)], Malignant=1-barOut$Nonmalignant, Clone.2=(barOut$clone.2+barOut$clone.3+barOut$clone.4+barOut$clone.5+barOut$clone.6), Clone.3=(barOut$clone.3+barOut$clone.4+barOut$clone.5))
colnames(barOutIn)[1:3]=c('Clone.4', 'Clone.5', 'Clone.6')
barOutIn9495=barOutIn[,c(4,5,6,1,2,3)]
barOutIn9495 = barOutIn9495[which(rownames(barOutIn9495) %in% as.numeric(gsub('X', '', colnames(rnaScale9495)[-1]))),]
rnaScale9495 = rnaScale9495[,-ncol(rnaScale9495)]
# }}}
# normalized SF10711 expression matrix - rnaScale10711
# {{{
rna=fread('/mnt/bdata/@patrick/SF10711/rna.seq/expression_matrices/Normalized_read_counts_using_RUVg_ERCC_K10Factors.csv', drop=seq(2,6))
colnames(rna)[1]='Gene'
meanExpr=apply(rna[,-seq(1,6)], 1, var)
rna=rna[order(meanExpr, decreasing=T),]
rnarSums=apply(rna[,-1], 1, sum)
rna=rna[-which(rnarSums<quantile(rnarSums,0.15)),]
rna=rna[-grep('^ERCC-', rna$Gene),]
rnarSumsN=apply(rna[,-1], 1, sum)
rna=as.data.frame(rna)
rnaScale10711 = rna[,-1]
# correct line below and make sure that the sections line up as expected
colnames(rnaScale10711)=gsub('SF10711)', '', colnames(rnaScale10711))
rnaScale10711=data.frame(Gene=rna$Gene, rnaScale10711)
# }}}
# clone 1 abundancd SF10711 - barOutIn10711
# {{{
barOut=barOutIn=read.csv('~/@patrick/SF10711/figures/tables/sf10711_cellular_abundance.csv')
rnaScale10711 = rnaScale10711[, c(1, which(colnames(rnaScale10711) %in% barOut$ampseq.sample))]
barOutIn10711=data.frame(barOut[,-c(1,2,3,4)], Malignant=1-barOut$nonmalignant, Clone.2=(barOut$clone.2+barOut$clone.3))
# }}}
# correlation analysis
# {{{
cor9495 = apply(rnaScale9495[,-1], 1, function(x) cor(x, barOutIn9495$Malignant))
names(cor9495) = rnaScale9495$Gene
cor10711 = apply(rnaScale10711[,-1], 1, function(x) cor(x, barOutIn10711$Malignant))
names(cor10711) = rnaScale10711$Gene
inter = (intersect(names(cor9495), names(cor10711)))
# }}}
# calculate p-values from correlation
# {{{
cor10711pVal = apply(rnaScale10711[,-1], 1, function(x) cor.test(x, barOutIn10711$Malignant)[['p.value']])
cor9495pVal = apply(rnaScale9495[,-1], 1, function(x) cor.test(x, barOutIn9495$Malignant)[['p.value']])
# }}}
# generate output table
# {{{
pValAll = data.frame(
    gene = inter,
    cor9495 = cor9495[match(inter, names(cor9495))],
    cor10711 = cor10711[match( inter, names(cor10711))],
    p9495 = cor9495pVal[match(inter, names(cor9495))], 
    p10711 = cor10711pVal[match(inter, names(cor10711))]    
)
pValAll = cbind(
    pValAll,
    cor9495Percentile =  rank(abs(pValAll$cor9495))/nrow(pValAll),
    cor10711Percentile = rank(abs(pValAll$cor10711))/nrow(pValAll)
)
# }}}
# combining p-values using fisher's method
# {{{
pValAll = cbind(pValAll, joint_pval = unlist(apply(pValAll[,c(2,3)], 1, function(x) fisher.transform(unlist(x)))))
pValAll = cbind(pValAll, joint_cor = unlist(apply(pValAll[,c(2,3)], 1, function(x) fisher.transformtocor(unlist(x)))))
signVec = rep(NA, nrow(pValAll))
signVec[pValAll$cor9495 > 0 & pValAll$cor10711 > 0] = 'positive'
signVec[pValAll$cor9495 < 0 & pValAll$cor10711 < 0] = 'negative'
pValAll$joint_qval = qvalue::qvalue(pValAll$joint_pval)$qvalues
pValAll = cbind(pValAll, 
                signVec, 
                sigVecBC = pValAll$joint_pval < (.05 / (nrow(pValAll) * 100)), 
                sigVecFDR = pValAll$joint_qval < .01
)
write.csv(pValAll, file = paste0('integration_figure_correlations_', cor_method, '_', method, '.csv'), row.names = F)
write.csv(pValAll, file = 'integration_figure_correlations__correlation_fisher.csv', row.names = F)
head(pValAll[order(pValAll$joint_pval, decreasing=F),], 100)
head(pValAll[order(pValAll$ cor10711Percentile, decreasing=T),], 100)
head(pValAll[order(pValAll$ cor10711Percentile, decreasing=T),], 100)
thresh = 0.3
thresh2 = 0.3
geneL = list(sigPosBC = pValAll[which(pValAll$signVec == 'positive' & pValAll$sigVecBC & pValAll$cor9495Percentile > thresh & pValAll$cor10711Percentile > thresh2),]$gene,
             sigNegBC = pValAll[which(pValAll$signVec == 'negative' & pValAll$sigVecBC & pValAll$cor9495Percentile > thresh & pValAll$cor10711Percentile > thresh2),]$gene,
             sigPosFDR = pValAll[which(pValAll$signVec == 'positive' & pValAll$sigVecFDR & pValAll$cor9495Percentile > thresh & pValAll$cor10711Percentile > thresh2),]$gene,
             sigNegFDR = pValAll[which(pValAll$signVec == 'negative' & pValAll$sigVecFDR & pValAll$cor9495Percentile > thresh & pValAll$cor10711Percentile > thresh2),]$gene
             )
enrichClust = enrichment_man(geneL, pValAll$gene, '/home/shared/genesets/genesets_slim')
head(enrichClust[order(enrichClust$sigPosFDR, decreasing = FALSE),])
head(enrichClust[order(enrichClust$sigNegFDR, decreasing = FALSE),])
# }}}
# import correlation and enrichment dataframes
# {{{
pValAll = data.frame(fread('~/code/git/integration_SF10711_SF9495/integration_figure_correlations_correlation_fisher.csv'))
cor.test(pValAll$cor9495, pValAll$cor10711)
2*pt(62.458, 15286, lower=F, log=T)
sets = c('MOSET6798', 'MOSET6718', 'MOSET6814', 'MOSET6760',
         'MOSET24', 'MOSET25', 'MOSET3', 'MOSET151')
geneL = list(sigPosBC = pValAll[which(pValAll$signVec == 'positive' & pValAll$sigVecBC & pValAll$cor9495>0.3 & pValAll$cor10711>0.3]$gene,
             sigNegBC = pValAll[which(pValAll$signVec == 'negative' & pValAll$sigVecBC & pValAll$cor9495<(-0.3) & pValAll$cor10711<(-0.3)),]$gene,
             sigPosFDR = pValAll[which(pValAll$signVec == 'positive' & pValAll$sigVecFDR & pValAll$cor9495>0.3 & pValAll$cor10711>0.3),]$gene,
             sigNegFDR = pValAll[which(pValAll$signVec == 'negative' & pValAll$sigVecFDR & pValAll$cor9495<(-0.3) & pValAll$cor10711<(-0.3)),]$gene
             )
sigI=which(colnames(pValAll) == 'sigVecBC')
enrichClust = enrichment_man(geneL, pValAll$gene, '/home/shared/genesets/genesets_slim')
write.csv(enrichClust, file = 'integration_figure_enrichments_our_sets_correlation_fisher_cutoff_0.3rel.csv', row.names = F)
enrichClust = enrichment_man_broad(geneL, pValAll$gene, broadSets)
write.csv(enrichClust, file = 'integration_figure_enrichments_broad_sets_correlation_fisher_cutoff_0.3rel.csv', row.names = F)
enrichClust =data.frame(fread('integration_figure_enrichments_our_sets_correlation_fisher_cutoff_0.3rel.csv'))
enrichClust = enrichClust[, -grep('^sig.*BC$', colnames(enrichClust))]
protein_aliases = geneL[[1]]
# }}}
# colors for dotplot
# {{{
sigCol = rep('c2', nrow(pValAll))
sigCol[which(pValAll$signVec == 'positive' & pValAll[,sigI] & pValAll$sigVecBC & pValAll$cor9495Percentile > thresh & pValAll$cor10711Percentile > thresh2)] = 'c1'
sigCol[which(pValAll$signVec == 'negative' & pValAll[,sigI] & pValAll$sigVecBC & pValAll$cor9495Percentile > thresh & pValAll$cor10711Percentile > thresh2)] = 'c3'
genesEnrich = list(pos = pValAll$gene[which(pValAll$signVec == 'positive' & pValAll[,sigI] & pValAll$sigVecBC & pValAll$cor9495>0.3 & pValAll$cor10711>0.3)], 
                   neg = pValAll$gene[which(pValAll$signVec == 'negative' & pValAll[,sigI] & pValAll$sigVecBC & pValAll$cor9495<(-0.3) & pValAll$cor10711<(-0.3))])
setNames = enrichClust$SetName[match(sets, enrichClust$SetID)]
# }}}
# set up enrichment plot
# {{{
enrichPlot = reshape2::melt(enrichClust[which(enrichClust$SetID %in% sets), c(2, seq(8,ncol(enrichClust)))], id.var = 'SetName')
enrichPlot$SetName[1] = 'Verhaak: classical subtype'
enrichPlot$SetName[2] = 'Ashkan: E14 radial glia'
enrichPlot$SetName[3] = 'Engler: inflitrating monocytes (astro)'
enrichPlot$SetName[4] = 'Fietz: extracellular matrix'
enrichPlot$SetName[5] = 'G2C: clathrin-coated vesicles'
enrichPlot$SetName[6] = 'G2C: synaptosome'
enrichPlot$SetName[7] = 'ABA: neurons'
enrichPlot$SetName[8] = 'Stanford: cerebral cortex'
enrichPlot$SetName[9:16] = enrichPlot$SetName[1:8]
rownames(enrichPlot)=seq(1,nrow(enrichPlot))
enrichPlot = enrichPlot[c(seq(1,4), 14,15,16,13),]
enrichPlot$value = -log10(enrichPlot$value)
enrichPlot$SetName = factor(enrichPlot$SetName, levels = enrichPlot$SetName[seq(length(setNames), 1)], ordered=T)
# }}}
# ggplots
# {{{
shapeL = rep(19, nrow(pValAll))
sizeL = rep(1, nrow(pValAll))
strokeL = rep(1, nrow(pValAll))
shapeL[grep('AKR1C3', pValAll$gene)]=8
sizeL[grep('AKR1C3', pValAll$gene)]=5
strokeL[grep('AKR1C3', pValAll$gene)]=3
dotenrich = ggplot(enrichPlot, aes(x=SetName, y=value, fill=variable)) +
    geom_bar(stat = "identity")+
    coord_flip() +
    scale_fill_manual(values = c('#1f78b4',  '#e31a1c')[seq(2,1)])+
    geom_text(aes(label = SetName, y=2), position = position_stack(vjust = 0.5), size = 13, hjust=0) +
    theme_classic() +scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    labs(y = '-Log10 q-value', x = '', title = 'Enrichments of highly\ncorrelated genes')+
    oldham_theme() +
    theme(legend.position = 'none', axis.text.y = element_blank(), axis.text.x =  element_text(size=30, color='black', family='NimbusSan', margin=margin(t=20)),plot.title = element_text(size=45,face="bold", hjust=.5, family='NimbusSan', margin=margin(t=-20, b=10)))
dotplotgg = ggplot(pValAll, aes(x = cor9495, y = cor10711, color = sigCol)) +
    geom_point(shape = shapeL, size = sizeL, stroke =strokeL) +
    theme_classic() +
    oldham_theme() +
    theme(legend.position = 'none') +
    labs(x = 'Correlation in case one', y = 'Correlation in case two', title = 'Correlation to Clone 1\nin both cases', subtitle = "Pearson's R = 0.45, p-value = 1E-300") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(method = "lm", se = FALSE, color = 'gray') +
    theme(legend.position = 'none', axis.text.x =  element_text(size=30, color='black', family='NimbusSan', margin=margin(t=20)), axis.text.y =  element_text(size=30, color='black', family='NimbusSan', margin=margin(t=20, r =20)),plot.title = element_text(size=45,face="bold", hjust=.5, family='NimbusSan', margin=margin(t=-20, b=10))) +
    scale_color_manual(values = c('#1f78b4', '#1a1a1a', '#e31a1c')[seq(3,1)])+
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 1), expand = c(0, 0))
pdf('integration_figure_final_correlation_fisher0.3_bc.pdf', height=13*.6, width=13*1.5)
grid.arrange(dotplotgg, dotenrich, ncol=2)
dev.off()
# }}}
# STRINGdb portion
# download PPI network from database
# {{{
file_version = '11.5'
species = 9606
score_threshold = 900
protocol = 'https'
stable_url = 'https://string-db.org/'
url <- paste(protocol, "://stringdb-static.org/download/protein.physical.links.v", file_version, "/", species, ".protein.physical.links.v", file_version, ".txt.gz", sep="")
temp = downloadAbsentFile(url)
PPI <- read.table(temp, sep = " ", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
PPIselected = PPI
if(length(score_threshold)!=0) PPIselected <- PPI[PPI$combined_score >= score_threshold,]
myg = graph.data.frame(PPIselected,FALSE)
graph = myg
# }}}
# remap gene names to STRING_id
# {{{
temp_df = data.frame(proteins=c(protein_aliases, 'MYC'))
aliases_type = 'takeFirst'
version = '11.5'
proteins = data.frame()
oD = getwd()
temp_df_mapped = map(temp_df, "proteins", removeUnmappedRows=TRUE, quiet=FALSE)
example1_mapped = (temp_df_mapped$STRING_id)
# }}}
# get PPI from chosen genes
# {{{
test =induced.subgraph(graph, which(V(graph)$name %in% temp_df_mapped$STRING_id))
clustTest = clusters(test)
clustChosen = names(table(clustTest$membership)[table(clustTest$membership)>5])
chosenGenes = names(clustTest$membership[clustTest$membership %in% clustChosen])
genesPClust = lapply(clustChosen, function(x) names(clustTest$membership[clustTest$membership %in% x]))
colVec = paste0('col', seq_along(genesPClust))
genesPcol = data.frame(genes = unlist(genesPClust), colors = unlist(lapply(seq_along(genesPClust), function(i) rep(colVec[i], length(genesPClust[[i]])))))
test = induced.subgraph(graph, which(V(graph)$name %in% temp_df_mapped$STRING_id[temp_df_mapped$STRING_id %in% chosenGenes]))
inter = ggnetwork(asNetwork(test), cell.jitter = .7, layout='fruchtermanreingold', repulse.rad = 400000) # possible: circle, kamadakawai, random
vNames = temp_df_mapped$proteins[match(inter$vertex.names, temp_df_mapped$STRING_id)]
vNamesClust = lapply(genesPClust, function(x) temp_df_mapped$protein[match(x, temp_df_mapped$STRING_id)])
keepers = c(
            'WNT3', 'COL1A1', 'SOX9', 'WNT5A','CDK2', 'VCAM1', 'POL2RG', 'SOX2', 'COL4A5' , 'SOX6', 'CDH7', 'CD44', 'CDKN2A', 'SOX6', 'MYC', 'CTNNB1',
            'CD72', 'BCL6',
            'SNRNP70', 'SNRPG', 'SNRPA',
            'ATM', 'POLG2', 'GATA1','XPC', 'CETN2', 'RBBP8',
            'NXF1', 'NUP85',
            'RING1', 'PHC1',
            'CENPA', 'NDC80',
            'RAD9A', 'RFC4'
    )
vNames[!(vNames %in% keepers)]=''
genesPcol = genesPcol[match(inter$vertex.names, genesPcol$genes),]
# }}}
# PPI enrichments
# {{{
enrichClust = enrichment_man(vNames, pValAll$gene, '/home/shared/genesets/genesets_slim')
write.csv(enrichClust, file = 'integration_figure_stringenrichments_our_sets_correlation_fisher_cutoff_0.3rel.csv', row.names = F)
enrichClust = enrichment_man_broad(vNames, pValAll$gene, broadSets)
write.csv(enrichClust, file = 'integration_figure_stringenrichments_broad_sets_correlation_fisher_cutoff_0.3rel.csv', row.names = F)
# }}}
# read in enrichments
enrichClustBroad = data.frame(fread('integration_figure_stringenrichments_broad_sets_correlation_fisher_cutoff_0.3rel.csv'))
# select genesets
# {{{
sets = c('M7847', 'M139', # cluster 1
         'M25144', 'M12336', # cluster 2
        'M13559', 'M9529',  # cluster 3
        'M11563', 'M13636', # cluster 4
        'M725', 'M17706', # cluster 5
        'M14738', 'M17402', # cluster 6
        'M25595', 'M17499', # cluster 7
        'M11153', 'M19381' # cluster 8
        )
setNames = enrichClustBroad$SetName[match(sets, enrichClustBroad$SetID)]
enrichPlot = reshape2::melt(enrichClustBroad[which(enrichClustBroad$SetID %in% sets), c(2, seq(8,ncol(enrichClustBroad)))], id.var = 'SetName')
enrichPlot$SetName = factor(enrichPlot$SetName, levels = setNames, ordered=T)
enrichPlot = enrichPlot[order(enrichPlot$SetName),]
rownames(enrichPlot)=seq(1,nrow(enrichPlot))
cV =c()
delta= 8
j = 1
for(i in seq_len(delta)){
    cV = c(cV, j, j + delta)
    j = j + delta*2 +1
}
enrichPlot = enrichPlot[cV, ]
enrichPlot$SetName = factor(enrichPlot$SetName, levels = setNames[seq(length(setNames), 1)], ordered=T)
enrichPlot$value = -log10(enrichPlot$value)
enrichPlot$SetName = capitalize(tolower(gsub('_', ' ', enrichPlot$SetName)))
enrichPlot$SetName = c('Reactome: signaling by WNT',
                       'PID: MYC pathway',
                       'GOBP: mononuclear cell differentiation',
                       'GOBP: regulation of immune response',
                       'GOBP: RNA splicing',
                       'GOBP: RNA processing',
                       'Kauffmann: DNA repair',
                       'GOBP: response to DNA damage',
                       'Reactome: transport of mRNA to cytplasm',
                       'GOCC: nuclear pore',
                       'GOCC: nuclear ubiquitin ligase',
                       'GOCC: PRC1 comples',
                       'GOCC: nuclear kinetichore',
                       'GOCC: centromeric region',
                       'Reactome: activation of ATR',
                       'Reactome: G2/M checkpoints'
                       )
enrichPlot$SetName = factor(enrichPlot$SetName, levels = enrichPlot$SetName[seq(length(setNames), 1)], ordered=T)
# }}}
# ploting STRINGdb network and enirchments
# {{{
labelCols = rep(brewer.pal(n=8, 'Set1'), times=unlist(lapply(vNamesClust,  length)))
names(labelCols) = unlist(vNamesClust)
labelCols = labelCols[match(vNames, names(labelCols))]
labelCols[is.na(labelCols)]='black'
labelCols = as.character(labelCols)
stringgg = ggplot(inter, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(), color = "black", size = 2) +
    geom_nodes(aes(color = genesPcol$colors), size=5)+
    geom_nodelabel_repel(aes(label = vNames, color=vNames),
                        color = 'black',
                        fill='white',
                        size=3,
                        segment.size=0,
                        family = 'NimbusSan',
                        label.padding = unit(0, 'lines'),
                        point.padding= NA,
                        max.overlaps=10) +
    scale_color_manual(values = c('#4daf4a', '#984ea3', '#ff7f00', '#d3d300', '#a65628', '#f781bf', '#999999', '#00d3c7'))+
    theme_blank() +
    oldham_theme() +
    labs(x ='', y='', title = 'STRING PPI of highly correlated genes', subtitle = paste0(length(unlist(vNamesClust))-1, ' proteins, ', nrow(inter), ' interactions (P < 2E-16)')) +
    theme(legend.position = 'none', axis.text.y = element_blank(), axis.text.x = element_blank(),plot.title = element_text(size=45,face="bold", hjust=.5, family='NimbusSan', margin=margin(t=-20, b=10)), axis.line.x = element_blank(), axis.line.y = element_blank())

enrichgg = ggplot(enrichPlot, aes(x=SetName, y=value, fill=variable)) +
    geom_bar(stat = "identity")+
    coord_flip() +
    scale_fill_manual(values = c('#4daf4a', '#984ea3', '#ff7f00', '#d3d300', '#a65628', '#f781bf', '#999999', '#00d3c7'))+
    geom_text(aes(label = SetName, y=2), position = position_stack(vjust = 0.5), size = 13, hjust=0) +
    theme_classic() +scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    labs(y = '-Log10 q-value', x = '', title = 'Enrichments of STRING PPI clusters')+
    oldham_theme() +
    theme(legend.position = 'none', axis.text.y = element_blank(), axis.text.x =  element_text(size=30, color='black', family='NimbusSan', margin=margin(t=20)),plot.title = element_text(size=45,face="bold", hjust=.5, family='NimbusSan', margin=margin(t=-20, b=10)))
pdf('integration_figure_final_string_correlation_fisher0.3_bc.pdf', height=13*1.2, width=30*1)
grid.arrange(stringgg, enrichgg, ncol=2)
dev.off()
# }}}
