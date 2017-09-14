library(ggplot2)
library(dplyr)
library(magrittr)
library(markerGeneProfile)
devtools::load_all()
set.seed(1)

StanleyC = data.frame(Gene.Symbol = rn(StanleyCExp),Probe= paste('p',1:nrow(StanleyCExp)), StanleyCExp)
StanleyC = data.frame(Gene.Symbol = rn(StanleyCExp), StanleyCExp)


genes = neuroExpressoAnalysis::mouseMarkerGenesCombined$Cortex['Oligo']

groups = StanleyCMeta$Profile
StanleyCEstimate = mgpEstimate(StanleyC,
                                    genes= genes,
                                    geneColName="Gene.Symbol",
                                    groups = groups,
                                    outlierSampleRemove= FALSE,
                                    removeMinority=TRUE,
                                    seekConsensus = FALSE)

# rawUranova = rawUranova[match(StanleyCMeta$`Stanley ID`,rawUranova$StanleyID),]

frame = data.frame(oligoEstim = StanleyCEstimate$estimates$Oligo %>% scale01, 
                   groups = StanleyCEstimate$groups$Oligo %>% factor(levels=c('Control', 
                                                                              'Schizophrenia',
                                                                              'Bipolar',
                                                                              'Depression'))# ,
                   # UranovaCounts = rawUranova$GM
)

frame$ids = StanleyCMeta$`Stanley ID`
plot(frame$oligoEstim,frame$UranovaCounts)
# cor(frame$oligoEstim, frame$UranovaCounts,method='spearman')

meanFrame = frame %>% group_by(groups) %>% summarise(mean(oligoEstim)) %>% 
    mutate(x1 = as.numeric(groups)-0.3, x2 = as.numeric(groups) +0.3)


wilcoxResults = sapply(c('Schizophrenia',
                   'Bipolar',
                   'Depression'), function(x){
                       a1 = frame$oligoEstim[frame$groups %in% 'Control']
                       a2 = frame$oligoEstim[frame$groups %in% x]
                       test = wilcox.test(a1,a2)
                       pValue = test$p.value
                       
                       W = unname(test$statistic)
                       
                       meanControl = mean(a1)
                       meanGroup = mean(a2)
                       
                       sdControl = sd(a1)
                       sdGroup = sd(a2)
                       return(c(pValue = pValue,W = W ,meanControl = meanControl, meanGroup = meanGroup,sdControl = sdControl, sdGroup = sdGroup ))
                   }) %>% t


pValues = wilcoxResults[,'pValue']

wilcoxResults %<>% round(digits = 3)

dir.create('analysis//03.MarkerGeneProfiles/tables', showWarnings=FALSE)
write.table(wilcoxResults, file =  'analysis//03.MarkerGeneProfiles/tables/stanleyCCounts.tsv',quote =FALSE,sep = '\t')

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
          strip.text.x = element_text(size = 19),
          legend.position='none') +
    geom_text(data=signifFrame,aes(label = text), size=10) + 
    scale_shape_manual(values = c(18, 2, 1, 0)) + 
    coord_cartesian(ylim = c(-0.03, 1.10))  + xlab('') + ylab('MGP estimation')

ggsave(filename='analysis//03.MarkerGeneProfiles/publishPlot/stanleyC.png',p,width=4.5,height=5,units='in')
ggsave(filename='analysis//03.MarkerGeneProfiles/publishPlot/stanleyC.svg',p,width=4.5,height=5,units='in')
ggsave(filename='analysis//03.MarkerGeneProfiles/publishPlot/stanleyC.pdf',p,width=4.5,height=5,units='in')

data = read.csv('data-raw/PlotExtract/plot.csv',header = F)
data %<>% filter(!V1<0.5) %>%
    arrange(V1) %>%
    mutate(groups = kmeans(V1,centers = 4)$cluster)

groupDictionary =  c('Control','Schizophrenia','Bipolar','Depression')
names(groupDictionary) = data$groups %>% unique

data %<>% mutate(groups = replaceElement(groups,groupDictionary)$newVect %>% factor(levels =  c('Control','Schizophrenia','Bipolar','Depression')))

# wilcoxResults = sapply(c('Schizophrenia',
#                          'Bipolar',
#                          'Depression'), function(x){
#                              a1 = data$V2[frame$groups %in% 'Control']
#                              a2 = data$V2[frame$groups %in% x]
#                              test = t.test(a1,a2)
#                              pValue = test$p.value
#                              
#                              W = unname(test$statistic)
#                              
#                              meanControl = mean(a1)
#                              meanGroup = mean(a2)
#                              
#                              sdControl = sd(a1)
#                              sdGroup = sd(a2)
#                              return(c(pValue = pValue,W = W ,meanControl = meanControl, meanGroup = meanGroup,sdControl = sdControl, sdGroup = sdGroup ))
#                          }) %>% t
# write.table(wilcoxResults, file =  'analysis//03.MarkerGeneProfiles/tables/stanleyCCounts.tsv',quote =FALSE,sep = '\t')

meanFrame = 
    data  %>% group_by(groups) %>% summarise(mean(V2)) %>% 
    mutate(x1 = as.numeric(groups)-0.3, x2 = as.numeric(groups) +0.3)

signifFrame = data.frame(groups = c('Schizophrenia',
                                    'Bipolar',
                                    'Depression'),
                         text = c('*','*','*'),
                         V2 = 110)
p = data %>% ggplot(aes(x = groups, y = V2, shape = groups)) +theme_cowplot(17) + 
    #geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    #geom_boxplot(width=0.1,fill = 'lightblue',outlier.size=0) + 
    geom_point(size = 3,position=position_jitter(width = 0.3), stroke  = 1.4)  + 
    #geom_jitter(size = 3)+ 
    # geom_point()+
    geom_segment(data=meanFrame, aes(x = x1, xend = x2, y = `mean(V2)`, yend = `mean(V2)`)) + 
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 19),
          axis.title.y = element_text(size = 14),
          legend.position='none') +
    geom_text(data=signifFrame,aes(label = text), size=10) + 
    scale_shape_manual(values = c(18, 2, 1, 0)) + 
    #coord_cartesian(ylim = c(-0.03, 1.10))  + 
    xlab('') + ylab(bquote('No of Oligodendroglial cells / 0.001 mm'^3))
ggsave(filename='analysis//03.MarkerGeneProfiles/publishPlot/stanleyCCounts.png',p,width=4.5,height=5,units='in')
ggsave(filename='analysis//03.MarkerGeneProfiles/publishPlot/stanleyCCounts.svg',p,width=4.5,height=5,units='in')
ggsave(filename='analysis//03.MarkerGeneProfiles/publishPlot/stanleyCCounts.pdf',p,width=4.5,height=5,units='in')

