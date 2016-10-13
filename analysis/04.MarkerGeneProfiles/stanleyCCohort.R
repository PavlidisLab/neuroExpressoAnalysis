devtools::load_all()
library(ggplot2)

StanleyC = data.frame(Gene.Symbol = rn(StanleyCExp), StanleyCExp)

genes = mouseMarkerGenes$Cortex['Oligo']
groups = StanleyCMeta$Profile
StanleyCEstimate = cellTypeEstimate(StanleyC,
                                    genes= genes,
                                    geneColName="Gene.Symbol",
                                    groups = groups,
                                    outlierSampleRemove= FALSE,
                                    removeNegatives=TRUE,
                                    seekConsensus = FALSE)

frame = data.frame(oligoEstim = StanleyCEstimate$estimates$Oligo %>% scale01, groups = StanleyCEstimate$groups$Oligo %>% factor(levels=c('Control', 
                                                                                                                  'Schizophrenia',
                                                                                                                  'Bipolar',
                                                                                                                  'Depression')))

meanFrame = frame %>% group_by(groups) %>% summarise(mean(oligoEstim)) %>% 
    mutate(x1 = as.numeric(groups)-0.3, x2 = as.numeric(groups) +0.3)


pValues = sapply(c('Schizophrenia',
                   'Bipolar',
                   'Depression'), function(x){
                       a1 = frame$oligoEstim[frame$groups %in% 'Control']
                       a2 = frame$oligoEstim[frame$groups %in% x]
                       wilcox.test(a1,a2)$p.value
                   })

marks = pValues %>% signifMarker

signifFrame = data.frame(groups = c('Schizophrenia',
                                    'Bipolar',
                                    'Depression'),
                         text = marks,
                         oligoEstim = 1.05)

p = frame %>% ggplot(aes(x = groups, y = oligoEstim, shape = groups)) +theme_cowplot(17) + 
    #geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    #geom_boxplot(width=0.1,fill = 'lightblue',outlier.size=0) + 
    geom_point(size = 3,position=position_jitter(width = 0.3), stroke  = 1.4)  + 
    #geom_jitter(size = 3)+ 
    # geom_point()+
    geom_segment(data=meanFrame, aes(x = x1, xend = x2, y = `mean(oligoEstim)`, yend = `mean(oligoEstim)`)) + 
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 14)) +
    geom_text(data=signifFrame,aes(label = text), size=10) + 
    scale_shape_manual(values = c(18, 2, 1, 0)) + 
    coord_cartesian(ylim = c(-0.03, 1.10))  + xlab('') + ylab('MGP estimation')

ggsave(filename='analysis//04.MarkerGeneProfiles/publishPlot/stanleyC.png',p,width=6,height=5,units='in')
