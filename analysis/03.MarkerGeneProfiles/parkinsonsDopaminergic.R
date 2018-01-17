# comparing dopaminergic cell proportions in substantia nigra samples
library(ogbox)
library(limma)
library(dplyr)
library(magrittr)
library(cowplot)
library(VennDiagram)
library(forcats)
library(purrr)
# library(markerGeneProfile)
devtools::load_all()

genes = mouseMarkerGenesCombined$Midbrain[
    !grepl('Microglia_',names(mouseMarkerGenesCombined$Midbrain))]
MoranDes = MoranParkinsonsMeta
# MoranDes %<>% mutate(patient = paste0(Disease, str_extract(Title,pattern='[0-9]+')))
list[geneDatMoran, expDatMoran] = sepExpr(MoranParkinsonsExp)
rownames(expDatMoran) = geneDatMoran$Gene.Symbol
# remove samples from superior frontal gyrus
expDatMoran = expDatMoran[(!MoranDes$Region %in% 'Superior frontal gyrus')]
MoranDes = MoranDes[(!MoranDes$Region %in% 'Superior frontal gyrus') ,]

ZhangDes = ZhangParkinsonsMeta %>% filter(brainRegion == 'Whole substantia nigra from postmortem brain')
list[geneDatZhang, expDatZhang] = sepExpr(ZhangParkinsonsExp)
expDatZhang = expDatZhang[ZhangParkinsonsMeta$GSM %in% ZhangDes$GSM]





expDats = list(Lesnick = LesnickParkinsonsExp,
               'Moran Lateral' = cbind(geneDatMoran, expDatMoran[MoranDes$Region %in% 'Lateral substantia nigra']),
               'Moran Medial' = cbind(geneDatMoran, expDatMoran[MoranDes$Region %in% "Medial substantia nigra"]),
               Zhang = cbind(geneDatZhang, expDatZhang))

groups = list(Lesnick = LesnickParkinsonsMeta$parkinson %>%
                  replaceElement(c('TRUE' = 'PD', 'FALSE' = 'control')) %$% newVector,
              'Moran Lateral' = MoranDes$Disease[MoranDes$Region %in% 'Lateral substantia nigra'],
              'Moran Medial' = MoranDes$Disease[MoranDes$Region %in% "Medial substantia nigra"],
              Zhang = ZhangDes$diseaseState)



# estimation for all
estimations = lapply(1:len(expDats),function(i){
    print(i)
    estimations =  mgpEstimate(exprData=expDats[[i]],
                               genes=genes,
                               geneColName='Gene.Symbol',
                               outlierSampleRemove=F,
                               groups=groups[[i]],
                               removeMinority = T,
                               PC = 1)
    
    wilcoxResults = estimations$estimates %>% sapply(function(x){
        x %<>% scale01
        grp = unique(groups[[i]])
        test = wilcox.test(x[groups[[i]] %in% grp[1]],x[groups[[i]] %in% grp[2]])
        p = test$p.value
        w = unname(test$statistic)
        
        controlMean = x[groups[[i]] %in% 'control'] %>% mean
        controlSD = x[groups[[i]] %in% 'control'] %>% sd
        nControl = sum(groups[[i]] %in% 'control')
        groupMean =  x[groups[[i]] %in% grp[!grp %in% 'control']] %>% mean
        groupSD =  x[groups[[i]] %in% grp[!grp %in% 'control']] %>% sd
        nGroup = sum(groups[[i]] %in% grp[!grp %in% 'control'])
        
        return(c(p = p, w=w , controlMean = controlMean, controlSD = controlSD,nControl=nControl,groupMean = groupMean,nGroup= nGroup, groupSD = groupSD))
        # p = wilcox.test(x[groups[[i]] %in% grp[1]],x[groups[[i]] %in% grp[2]])$p.value
    }) %>% t
    
    pVals = wilcoxResults[,'p']
    return(list(estimations = estimations,pVals = pVals,wilcoxResults = wilcoxResults))
})
names(estimations) = names(expDats)

dir.create('analysis/03.MarkerGeneProfiles/estimates',showWarnings = FALSE,recursive = TRUE)
saveRDS(estimations,file = 'analysis/03.MarkerGeneProfiles/estimates/parkinsonsEstimate.rds')

statsTable = estimations %>% sapply(function(x){
    x$wilcoxResults['Dopaminergic',]
}) %>% t %>% round(digits = 3)

write.table(statsTable, file =  'analysis//03.MarkerGeneProfiles/tables/parkinson.tsv',quote =FALSE,sep = '\t')

# dopaminergic gene counts
plotNames = sapply(1:len(estimations), function(i){
    geneCount = estimations[[i]]$estimations$rotations$Dopaminergic %>% nrow
    paste0(names(estimations)[i],'\n(n genes = ',geneCount,')')
})

# frame for plotting 
frames = lapply(1:len(estimations), function(i){
    geneCount = estimations[[i]]$estimations$rotations$Dopaminergic %>% nrow
    name =  paste0(names(estimations)[i],'\n(n genes = ',geneCount,')')
    
    frame = data.frame(parkinsons = estimations[[i]]$estimations$groups$Dopaminergic,
                       estimate = scale01(estimations[[i]]$estimations$estimates$Dopaminergic),
                       name = name,stringsAsFactors = FALSE)
})

masterFrame = rbindMult(list = frames)

pVals = estimations %>%
    purrr::map('pVals') %>%
    purrr::map_dbl('Dopaminergic') %>% 
    ogbox::signifMarker()

signifFrame = data.frame(markers = pVals,
                         x = 1.5,
                         y = 1.0,
                         name =frames %>% 
                             map('name') %>% map_chr(unique))



pEstimate = masterFrame %>%  ggplot(aes( y = estimate, x = parkinsons)) + 
    #geom_point(position= 'jitter',size=3) +
    facet_grid(~name)  +
    theme_cowplot(17) + 
    geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    geom_boxplot(width=0.2,fill = 'lightblue') + 
    # geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 13)) +
    coord_cartesian(ylim = c(-0.10, 1.10)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10)+
    xlab('') +
    ylab('Dopaminergic MGP estimation')

ggsave(plot=pEstimate,
       filename='analysis//03.MarkerGeneProfiles/publishPlot/dopaminergicEstimation.png',width=7,height=4.5,units='in')


# paper gene correlations -------------------
c('CDC42', 'FGF13',
  "HSPB1","SNCA",
  "MKNK2", "TF",
  "AMPH", "BEX1",
  "JMJD6", "NSF",
  "SUB1", "SV2B",
  "SYT1", "SNAP25",
  "STMN2", "RGS4",
  "SNX10", "PRKAR2B",
  "NEFL", "MDH1",
  "CHGB",
  "NFASC") %>% sort -> paperGenes

# merge different moran regions so you won't have too many plots
expDats = list(Lesnick = LesnickParkinsonsExp,
               'Moran' = cbind(geneDatMoran, expDatMoran),
               Zhang = cbind(geneDatZhang, expDatZhang))

groups = list(Lesnick = LesnickParkinsonsMeta$parkinson %>%
                  replaceElement(c('TRUE' = 'PD', 'FALSE' = 'control')) %$% newVector,
              Moran = MoranDes$Disease,
              Zhang = ZhangDes$diseaseState %>% 
                  replaceElement(c(Control = 'control', 'Parkinsons disease' = 'PD')) %$% newVector)

corPlotFrames = lapply(1:len(expDats),function(i){
    print(i)
    dopaEstim =  mgpEstimate(exprData=expDats[[i]],
                             genes=genes,
                             geneColName='Gene.Symbol',
                             outlierSampleRemove=F,
                             groups=groups[[i]],
                             removeMinority = T,
                             PC = 1)$estimates$Dopaminergic %>% scale01
    
    list[,paperGeneExp] = expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,] %>% sepExpr
    PC = paperGeneExp %>% t %>% prcomp(scale=TRUE) %$% x[,1]
    cor = (cor(PC,dopaEstim) > 0) - (cor(PC,dopaEstim) < 0)
    
    data.frame(disease = groups[[i]],
               paperGenePC = PC*(cor),
               estimation = dopaEstim,
               source = names(expDats)[i])
})
names(corPlotFrames) = names(expDats)
masterCorPlot = rbindMult(list = corPlotFrames )

groupCorrelations = lapply(1:len(expDats), function(i){
    controlCor = cor(corPlotFrames[[i]] %>% filter(disease == 'control') %$% estimation,
                     corPlotFrames[[i]] %>% filter(disease == 'control') %$% paperGenePC,
                     method='spearman') %>% abs
    
    PDCor = cor(corPlotFrames[[i]] %>% filter(disease == 'PD') %$% estimation,
                corPlotFrames[[i]] %>% filter(disease == 'PD') %$% paperGenePC,
                method='spearman') %>% abs
    
    c(control = controlCor, PD = PDCor )
})

names(groupCorrelations) = names(expDats)

annotation = data.frame(source = groupCorrelations %>% names,
                        label = c(paste0('Spearman\'s ρ:\n',
                                         'controls: ',format(groupCorrelations %>% map_dbl('control'),digits = 3),'\n',
                                         'PD: ', format(groupCorrelations %>% map_dbl('PD'),digits=3))))


onlySigFigures <- function(){
    # return a function responpsible for formatting the 
    # axis labels with a given number of decimals 
    function(x) as.char(x)
}

genePCestimationPlot = masterCorPlot %>% ggplot(aes(x =estimation , y = paperGenePC, color = disease)) + 
    facet_grid(~source) +
    scale_x_continuous(breaks= c(0,0.5,1),labels = onlySigFigures()) + 
    theme_cowplot(17) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 17))+
    geom_point(size=3) + 
    stat_smooth(method=lm, se=FALSE) + 
    scale_color_manual(values=c(viridis(5)[1],viridis(5)[4])) +
    xlab('Dopaminergic MGP estimation') +
    ylab('PC1 of PD signature gene expression') + 
    geom_text(data = annotation,
              aes(label = label,
                  y = min(masterCorPlot$paperGenePC),
                  x = max(masterCorPlot$estimation)),
              color = 'black',
              vjust= 0,
              hjust= 1,
              size = 5)

ggsave(filename = 'analysis//03.MarkerGeneProfiles/publishPlot/genePCestimation.png',
       plot= genePCestimationPlot,
       width=8.5,height=4.3,units='in')


# all genes correlation
allGeneCors = lapply(1:len(expDats), function(i){
    controlCor = cor(corPlotFrames[[i]]$estimation[corPlotFrames[[i]]$disease %in% 'control' ],
                     expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,] %>% 
                         sepExpr %>% {.[[2]][,corPlotFrames[[i]]$disease %in% 'control' ]} %>% t) %>% t
    
    controlCor = data.frame(Correlation = controlCor, disease = 'control', 
                            gene = expDats[[i]]$Gene.Symbol[expDats[[i]]$Gene.Symbol %in% paperGenes])
    
    PDcor =  cor(corPlotFrames[[i]]$estimation[corPlotFrames[[i]]$disease %in% 'PD' ],
                 expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,] %>% 
                     sepExpr %>% {.[[2]][,corPlotFrames[[i]]$disease %in% 'PD' ]} %>% t) %>% t
    
    PDcor = data.frame(Correlation = PDcor, disease = 'PD', 
                       gene = expDats[[i]]$Gene.Symbol[expDats[[i]]$Gene.Symbol %in% paperGenes])
    
    data.frame(rbind(PDcor,controlCor),source = names(expDats)[i])
})
allGeneCors %<>% rbindMult(list = .)

allGeneCors$disease %<>% factor(levels=c('control','PD'))

allGeneCors$gene %<>% fct_reorder(allGeneCors$Correlation,.desc=TRUE)


geneAllestimation = allGeneCors %>% ggplot(aes(x = gene, y = Correlation, color = disease)) + 
    facet_grid(~source) +
    geom_abline(slope = 0, intercept = 0 ,linetype =2 ) + 
    theme_cowplot(17) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 17))+
    geom_point(size=3, alpha = 0.8) + 
    scale_color_manual(values = c(viridis(5)[1],viridis(5)[4])) +
    xlab('') +
    ylab('Dopaminergic MGP-\nGene expression correlation') + 
    theme(axis.text.x  = element_text(angle= 90,vjust = 0.5, size = 13))

ggsave(filename = 'analysis//03.MarkerGeneProfiles/publishPlot/geneAllestimation.png',
       plot= geneAllestimation,
       width=13,height=4.3,units='in')


# Zhang other regions --------------------
# list[geneDatZhang, expDatZhang] = sepExpr(ZhangParkinsonsExp)
# 
# 
# expDats = list('Zhang SN' = cbind(geneDatZhang,expDatZhang[ZhangParkinsonsMeta$brainRegion %in% 'Whole substantia nigra from postmortem brain'] ),
#                'Zhang Putamen' =  cbind(geneDatZhang,expDatZhang[ZhangParkinsonsMeta$brainRegion %in% 'Putamen from postmortem brain'] ),
#                'Zhang Cortex' =  cbind(geneDatZhang,expDatZhang[ZhangParkinsonsMeta$brainRegion %in% 'Prefrontal cortex area 9'] ))
# 
# groups = list('Zhang SN' = ZhangParkinsonsMeta$diseaseState[ZhangParkinsonsMeta$brainRegion %in% 'Whole substantia nigra from postmortem brain'],
#               'Zhang Putamen' = ZhangParkinsonsMeta$diseaseState[ZhangParkinsonsMeta$brainRegion %in% 'Putamen from postmortem brain'],
#               'Zhang Cortex' =  ZhangParkinsonsMeta$diseaseState[ZhangParkinsonsMeta$brainRegion %in% 'Prefrontal cortex area 9'])
# 
# patients = list('Zhang SN' = ZhangParkinsonsMeta$patient[ZhangParkinsonsMeta$brainRegion %in% 'Whole substantia nigra from postmortem brain'],
#                 'Zhang Putamen' = ZhangParkinsonsMeta$patient[ZhangParkinsonsMeta$brainRegion %in% 'Putamen from postmortem brain'],
#                 'Zhang Cortex' =  ZhangParkinsonsMeta$patient[ZhangParkinsonsMeta$brainRegion %in% 'Prefrontal cortex area 9'])
# 
# estimations = lapply(1:len(expDats),function(i){
#     print(i)
#     estimations =  mgpEstimate(exprData=expDats[[i]],
#                                genes=genes,
#                                geneColName='Gene.Symbol',
#                                outlierSampleRemove=F,
#                                groups=groups[[i]],
#                                removeMinority = F,
#                                PC = 1)
#     
#     pVals = estimations$estimates %>% sapply(function(x){
#         grp = unique(groups[[i]])
#         p = wilcox.test(x[groups[[i]] %in% grp[1]],x[groups[[i]] %in% grp[2]])$p.value
#     })
#     return(list(estimations = estimations,pVals = pVals))
# })
# names(estimations) = names(expDats)
# # estimations %>% map('pVals') %>% map('Dopaminergic')
# 
# 
# plotNames = sapply(1:len(estimations), function(i){
#     geneCount = estimations[[i]]$estimations$rotations$Dopaminergic %>% nrow
#     paste0(names(estimations)[i],'\n(n genes = ',geneCount,')')
# })
# 
# # frame for plotting 
# frames = lapply(1:len(estimations), function(i){
#     geneCount = estimations[[i]]$estimations$rotations$Dopaminergic %>% nrow
#     name =  paste0(names(estimations)[i],'\n(n genes = ',geneCount,')')
#     
#     frame = data.frame(parkinsons = estimations[[i]]$estimations$groups$Dopaminergic,
#                        estimate = scale01(estimations[[i]]$estimations$estimates$Dopaminergic),
#                        name = name,stringsAsFactors = FALSE)
# })
# 
# masterFrame = rbindMult(list = frames)
# 
# pVals = estimations %>%
#     purrr::map('pVals') %>%
#     purrr::map_dbl('Dopaminergic') %>% 
#     ogbox::signifMarker()
# 
# signifFrame = data.frame(markers = pVals,
#                          x = 1.5,
#                          y = 1.0,
#                          name =frames %>% 
#                              map('name') %>% map_chr(unique))
# 
# 
# 
# pEstimate = masterFrame %>%  ggplot(aes( y = estimate, x = parkinsons)) + 
#     #geom_point(position= 'jitter',size=3) +
#     facet_grid(~name)  +
#     theme_cowplot(17) + 
#     geom_violin( color="#C4C4C4", fill="#C4C4C4") +
#     geom_boxplot(width=0.2,fill = 'lightblue') + 
#     # geom_point()+
#     theme(axis.text.x = element_text(angle=45, hjust = 1),
#           strip.text.x = element_text(size = 13)) +
#     coord_cartesian(ylim = c(-0.10, 1.10)) +
#     geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10)+
#     xlab('') +
#     ylab('Dopaminergic MGP estimation')
# 
# (pEstimate)
# 
# 
# corPlotFrames = lapply(1:len(expDats),function(i){
#     print(i)
#     dopaEstim =  mgpEstimate(exprData=expDats[[i]],
#                              genes=genes,
#                              geneColName='Gene.Symbol',
#                              outlierSampleRemove=F,
#                              groups=groups[[i]],
#                              removeMinority = F,
#                              PC = 1)$estimates$Dopaminergic %>% scale01
#     
#     list[,paperGeneExp] = expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,] %>% sepExpr
#     PC = paperGeneExp %>% t %>% prcomp(scale=TRUE) %$% x[,1] 
#     cor = (cor(PC,dopaEstim) > 0) - (cor(PC,dopaEstim) < 0)
#     
#     data.frame(disease = groups[[i]],
#                paperGenePC = PC*(cor),
#                estimation = dopaEstim,
#                source = names(expDats)[i],
#                patients = patients[[i]],
#                GSM = names(PC))
# })
# names(corPlotFrames) = names(expDats)
# 
# masterCorPlot = rbindMult(list = corPlotFrames )
# 
# groupCorrelations = lapply(1:len(expDats), function(i){
#     controlCor = cor(corPlotFrames[[i]] %>% filter(disease == 'control') %$% estimation,
#                      corPlotFrames[[i]] %>% filter(disease == 'control') %$% paperGenePC,
#                      method='spearman') %>% abs
#     
#     PDCor = cor(corPlotFrames[[i]] %>% filter(disease == 'PD') %$% estimation,
#                 corPlotFrames[[i]] %>% filter(disease == 'PD') %$% paperGenePC,
#                 method='spearman') %>% abs
#     
#     c(control = controlCor, PD = PDCor )
# })
# 
# names(groupCorrelations) = names(expDats)
# 
# annotation = data.frame(source = groupCorrelations %>% names,
#                         label = c(paste0('Spearman\'s ρ:\n',
#                                          'controls: ',format(groupCorrelations %>% map_dbl('control'),digits = 3),'\n',
#                                          'PD: ', format(groupCorrelations %>% map_dbl('PD'),digits=3))))
# 
# 
# genePCestimationPlot = masterCorPlot %>% ggplot(aes(x =estimation , y = paperGenePC, color = disease)) + 
#     facet_grid(~source) +
#     theme_cowplot(17) +
#     theme(strip.background = element_blank(),
#           strip.text = element_text(size = 17))+
#     geom_point(size=3) + 
#     stat_smooth(method=lm, se=FALSE) + 
#     scale_color_viridis(discrete=TRUE) +
#     xlab('Dopaminergic MGP estimation') +
#     ylab('PC1 of PD signature gene expression') + 
#     geom_text(data = annotation,
#               aes(label = label,
#                   y = min(masterCorPlot$paperGenePC),
#                   x = max(masterCorPlot$estimation)),
#               color = 'black',
#               vjust= 0,
#               hjust= 1,
#               size = 5)
# (genePCestimationPlot)
# 
# # all genes correlation
# allGeneCors = lapply(1:len(expDats), function(i){
#     controlCor = cor(corPlotFrames[[i]]$estimation[corPlotFrames[[i]]$disease %in% 'control' ],
#                      expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,] %>% 
#                          sepExpr %>% {.[[2]][,corPlotFrames[[i]]$disease %in% 'control' ]} %>% t) %>% t
#     
#     controlCor = data.frame(Correlation = controlCor, disease = 'control', 
#                             gene = expDats[[i]]$Gene.Symbol[expDats[[i]]$Gene.Symbol %in% paperGenes])
#     
#     PDcor =  cor(corPlotFrames[[i]]$estimation[corPlotFrames[[i]]$disease %in% 'PD' ],
#                  expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,] %>% 
#                      sepExpr %>% {.[[2]][,corPlotFrames[[i]]$disease %in% 'PD' ]} %>% t) %>% t
#     
#     PDcor = data.frame(Correlation = PDcor, disease = 'PD', 
#                        gene = expDats[[i]]$Gene.Symbol[expDats[[i]]$Gene.Symbol %in% paperGenes])
#     
#     data.frame(rbind(PDcor,controlCor),source = names(expDats)[i])
# })
# allGeneCors %<>% rbindMult(list = .)
# 
# allGeneCors$disease %<>% factor(levels=c('control','PD'))
# 
# allGeneCors$gene %<>% fct_reorder(abs(allGeneCors$Correlation),.desc=TRUE)
# 
# 
# geneAllestimation = allGeneCors %>% ggplot(aes(x = gene, y = Correlation, color = disease)) + 
#     facet_grid(~source) +
#     geom_abline(slope = 0, intercept = 0 ,linetype =2 ) + 
#     theme_cowplot(17) +
#     theme(strip.background = element_blank(),
#           strip.text = element_text(size = 17))+
#     geom_point(size=3, alpha = 0.8) + 
#     scale_color_viridis(discrete=TRUE) +
#     xlab('') +
#     ylab('Dopaminergic MGP-\nGene expression correlation') + 
#     theme(axis.text.x  = element_text(angle= 90,vjust = 0.5, size = 13))
# (geneAllestimation)
# 
# 
# 
# 
# # all genes correlation to SN -------------------------
# allGeneCors = lapply(1:len(expDats), function(i){
#     
#     snFrame = corPlotFrames[["Zhang SN"]][match(as.char(corPlotFrames[[i]]$patients ),
#                                                 as.char(corPlotFrames[["Zhang SN"]]$patients)) %>% trimNAs,]
#     thisFrame = corPlotFrames[[i]][match(as.char(corPlotFrames[['Zhang SN']]$patients ),
#                                          as.char(corPlotFrames[[i]]$patients)) %>% trimNAs,]
#     
#     controlCor = cor(snFrame$estimation[snFrame$disease %in% 'control' ],
#                      expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,as.char(thisFrame$GSM)] %>% 
#                          sepExpr %>% {.[[2]][,thisFrame$disease %in% 'control' ]} %>% t) %>% t
#     
#     controlCor = data.frame(Correlation = controlCor, disease = 'control', 
#                             gene = expDats[[i]]$Gene.Symbol[expDats[[i]]$Gene.Symbol %in% paperGenes])
#     
#     PDcor =  cor(corPlotFrames$`Zhang SN`$estimation[corPlotFrames[[i]]$disease %in% 'PD' ],
#                  expDats[[i]][expDats[[i]]$Gene.Symbol %in% paperGenes,] %>% 
#                      sepExpr %>% {.[[2]][,corPlotFrames[[i]]$disease %in% 'PD' ]} %>% t) %>% t
#     
#     PDcor = data.frame(Correlation = PDcor, disease = 'PD', 
#                        gene = expDats[[i]]$Gene.Symbol[expDats[[i]]$Gene.Symbol %in% paperGenes])
#     
#     data.frame(rbind(PDcor,controlCor),source = names(expDats)[i])
# })
# allGeneCors %<>% rbindMult(list = .)
# 
# allGeneCors$disease %<>% factor(levels=c('control','PD'))
# 
# allGeneCors$gene %<>% fct_reorder(abs(allGeneCors$Correlation),.desc=TRUE)
# 
# 
# geneAllestimation = allGeneCors %>% ggplot(aes(x = gene, y = Correlation, color = disease)) + 
#     facet_grid(~source) +
#     geom_abline(slope = 0, intercept = 0 ,linetype =2 ) + 
#     theme_cowplot(17) +
#     theme(strip.background = element_blank(),
#           strip.text = element_text(size = 17))+
#     geom_point(size=3, alpha = 0.8) + 
#     scale_color_viridis(discrete=TRUE) +
#     xlab('') +
#     ylab('Dopaminergic MGP-\nGene expression correlation') + 
#     theme(axis.text.x  = element_text(angle= 90,vjust = 0.5, size = 13))
# (geneAllestimation)
# 
# ####################
# masterCorPlot = rbindMult(list = corPlotFrames )
# 
# groupCorrelations = lapply(1:len(expDats), function(i){
#     controlCor = cor(corPlotFrames[[i]] %>% filter(disease == 'control') %$% estimation,
#                      corPlotFrames[[i]] %>% filter(disease == 'control') %$% paperGenePC,
#                      method='spearman') %>% abs
#     
#     PDCor = cor(corPlotFrames[[i]] %>% filter(disease == 'PD') %$% estimation,
#                 corPlotFrames[[i]] %>% filter(disease == 'PD') %$% paperGenePC,
#                 method='spearman') %>% abs
#     
#     c(control = controlCor, PD = PDCor )
# })
# 
# names(groupCorrelations) = names(expDats)
# 
# annotation = data.frame(source = groupCorrelations %>% names,
#                         label = c(paste0('Spearman\'s ρ:\n',
#                                          'controls: ',format(groupCorrelations %>% map_dbl('control'),digits = 3),'\n',
#                                          'PD: ', format(groupCorrelations %>% map_dbl('PD'),digits=3))))
# 
# 
# genePCestimationPlot = masterCorPlot %>% ggplot(aes(x =estimation , y = paperGenePC, color = disease)) + 
#     facet_grid(~source) +
#     theme_cowplot(17) +
#     theme(strip.background = element_blank(),
#           strip.text = element_text(size = 17))+
#     geom_point(size=3) + 
#     stat_smooth(method=lm, se=FALSE) + 
#     scale_color_viridis(discrete=TRUE) +
#     xlab('Dopaminergic MGP estimation') +
#     ylab('PC1 of PD signature gene expression') + 
#     geom_text(data = annotation,
#               aes(label = label,
#                   y = min(masterCorPlot$paperGenePC),
#                   x = max(masterCorPlot$estimation)),
#               color = 'black',
#               vjust= 0,
#               hjust= 1,
#               size = 5)
# (genePCestimationPlot)