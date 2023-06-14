# violin plot figure
# glasso
# need to load in Robj from previous workspace
# {{{
library('data.table')
library('ggplot2')
library('ggsignif')
library('kSamples')

base::load('~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/workspace_after_bootstrap.Robj')
base::load('~/@patrick/SF10711/integration_analysis/lasso_modeling/lasso/cumulative/workspace_after_bootstrapPerm.Robj')
clones=c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
tvalReal=tvalPerm=list()
for(j in seq_along(clones)){
    for(i in seq_along(coefOut)){
        x=unlist(lapply(coefOut[[i]], function(x) (names(x)==clones[j] & length(x)==1)[1]))
        if(i == 1) {
            tvalReal[[j]]=unlist(tvalOut[[i]][x])
            tvalPerm[[j]]=unlist(tvalOutPerm[[i]][x])
        } else {
            tvalReal[[j]]=c(tvalReal[[j]], unlist(tvalOut[[i]][x]))
            tvalPerm[[j]]=c(tvalPerm[[j]], unlist(tvalOutPerm[[i]][x]))
        }
    }
}
# }}}
# add names
# {{{
names(tvalReal) = clones
names(tvalPerm) = clones
tvalRealDf=reshape2::melt(tvalReal)
tvalPermDf=reshape2::melt(tvalPerm)
tvalRealDf = data.table(tvalRealDf, Permuted=rep('Real', nrow(tvalRealDf)))
tvalPermDf = data.table(tvalPermDf, Permuted=rep('Permuted', nrow(tvalPermDf)))
tvalPermDf$value[!(tvalPermDf$L1 == 'Clone.2')]=tvalPermDf$value[!(tvalPermDf$L1 == 'Clone.2')]
tvalPlot = rbind(tvalRealDf, tvalPermDf)
tvalPlot$L1 = gsub('\\.', ' ', tvalPlot$L1)
tvalPlot$Permuted = factor(tvalPlot$Permuted, levels=c('Real', 'Permuted'))
tvalPlot$L1 = gsub('Malignant', 'Clone 1', tvalPlot$L1)
tvalPlot$L1=factor(tvalPlot$L1, levels=c('Clone 1', 'Clone 2', 'Clone 3', 'Clone 4', 'Clone 5'))
# }}}
# do statistical testing
# {{{
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 2' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 2' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 2' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Permuted'][samp3k])
# }}}
# plot violin plots
# {{{
p1 = ggplot (tvalPlot, aes(x=L1, y=value, fill=Permuted)) + 
    geom_violin() +
    theme_classic() +
    oldham_theme() +
    scale_y_continuous(limits=c(-20, 20), breaks=seq(-20,20,20)) +
    labs(x='', y='T-values', fill='', title='T-values of lasso model') + # , subtitle='P-values from two-sided Anderson-Darling test') +
    scale_fill_manual(values=c('black', 'white')) +
    theme(legend.spacing.y = unit(1, 'cm'),
            axis.text.x = element_text(angle = 90)) +
    geom_signif(
        y_position = c(rep(18,5)),
        xmin = seq(0.8, 4.8,1),
        xmax = seq(1.2, 5.2,1),
        annotation = c('5E-92','3E-75', '3E-224', '2E-245', '5E-77'),
        tip_length = 0,
        size=2,
        textsize=10,
        vjust=-.6,
        family='NimbusSan'
    ) +
    guides(fill = guide_legend(byrow = TRUE))
# }}}
# lasso
# need to load in Robj from previous workspace
# {{{
base::load('~/@patrick/SF10711/integration_analysis/lasso_modeling/glasso/cumulative/workspace_after_bootstrap.Robj')
base::load('~/@patrick/SF10711/integration_analysis/lasso_modeling/glasso/cumulative/workspace_after_bootstrapPerm.Robj')
clones=c('Clone.1', 'Clone.2', 'Clone.3', 'Clone.4', 'Clone.5')
tvalReal=tvalPerm=list()
for(j in seq_along(clones)){
    for(i in seq_along(coefOut)){
        x=unlist(lapply(coefOut[[i]], function(x) (names(x)==clones[j] & length(x)==1)[1]))
        if(i == 1) {
            tvalReal[[j]]=unlist(tvalOut[[i]][x])
            tvalPerm[[j]]=unlist(tvalOutPerm[[i]][x])
        } else {
            tvalReal[[j]]=c(tvalReal[[j]], unlist(tvalOut[[i]][x]))
            tvalPerm[[j]]=c(tvalPerm[[j]], unlist(tvalOutPerm[[i]][x]))
        }
    }
}
# }}}
# add names
# {{{
names(tvalReal) = clones
names(tvalPerm) = clones
tvalRealDf=reshape2::melt(tvalReal)
tvalPermDf=reshape2::melt(tvalPerm)
tvalRealDf = data.table(tvalRealDf, Permuted=rep('Real', nrow(tvalRealDf)))
tvalPermDf = data.table(tvalPermDf, Permuted=rep('Permuted', nrow(tvalPermDf)))
tvalPlot = rbind(tvalRealDf, tvalPermDf)
tvalPlot$L1 = gsub('\\.', ' ', tvalPlot$L1)
tvalPlot$Permuted = factor(tvalPlot$Permuted, levels=c('Real', 'Permuted'))
tvalPlot$L1 = gsub('Malignant', 'Clone 1', tvalPlot$L1)
tvalPlot$L1=factor(tvalPlot$L1, levels=c('Clone 1', 'Clone 2', 'Clone 3', 'Clone 4', 'Clone 5'))
# }}}
# statistical testing
# {{{
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Permuted'][samp3k])$ad[2,3]
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 1' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 2' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 2' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 2' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 3' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 4' & tvalPlot$Permuted=='Permuted'][samp3k])
samp3k=sample(seq(1, length(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'])), 3000)
ad.test(tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Real'][samp3k], tvalPlot$value[tvalPlot$L1=='Clone 5' & tvalPlot$Permuted=='Permuted'][samp3k])
# }}}
# plot violin plots
# {{{
p2 = ggplot (tvalPlot, aes(x=L1, y=value, fill=Permuted)) + 
    geom_violin() +
    theme_classic() +
    oldham_theme() +
    scale_y_continuous(limits=c(-20, 20), breaks=seq(-20,20,20)) +
    labs(x='', y='T-values', fill='', title='T-values of group lasso model') +# , subtitle='P-values from two-sided Anderson-Darling test') +
    scale_fill_manual(values=c('black', 'white')) +
    theme(legend.spacing.y = unit(1, 'cm'),
            axis.text.x = element_text(angle = 90)) +
    geom_signif(
        y_position = c(rep(18,5)),
        xmin = seq(0.8, 4.8,1),
        xmax = seq(1.2, 5.2,1),
        annotation = c('9E-30','5E-45', '5E-108', '1E-165', '1E-21'),
        tip_length = 0,
        size=2,
        textsize=10,
        vjust=-.6,
        family='NimbusSan'
    ) +
    guides(fill = guide_legend(byrow = TRUE))

library('gridExtra')
p1 = p1 + 
    theme(axis.text.x=element_blank())
pdf('~/@patrick/SF10711/figures/tvalue_lasso.pdf', height=13*2, width=13)
grid.arrange(arrangeGrob(p1, p2, nrow=2), ncol=1)
dev.off()
p1_alt = p1 + theme(legend.position = 'none', plot.margin = margin(b = 0)) + scale_y_continuous(limits=c(-20, 20), breaks=seq(-20,20,20))
p2_alt = p2 + theme(legend.position = 'none',  plot.margin = margin(t = 0))+ scale_y_continuous(limits=c(-20, 20), breaks=seq(-20,20,20))
pdf('~/@patrick/SF10711/figures/tvalue_lasso_no_leg.pdf', height=13, width=13*.8)
grid.arrange(arrangeGrob(p1_alt, p2_alt, nrow=2, heights=c(4/10,5/10), padding = 0), ncol=1)
dev.off()
# }}}
