# devtools::load_all()
# # library(markerGeneProfile)
# 
# justMicroarray = pickMarkersAll('analysis/01.SelectGenes/FinalGenes1/PyramidalDeep/')
# withSingleCells = pickMarkersAll('analysis/01.SelectGenes/FinalGenes3/PyramidalDeep/')
# withSingleCells = pickMarkersAll('analysis/01.SelectGenes/FinalGenes3/PyramidalDeepNoSingleCellUnlessYouHaveTo//')
# 
# justSingleCells = pickMarkersAll('analysis/01.SelectGenes/QuickJustSingleCell/')
# 
# markerGenes = names(withSingleCells$All) %>% sapply(function(x){
#     c(withSingleCells$All[[x]],justSingleCells$PyramidalDeep[[x]]) %>% unique
# },simplify=FALSE)
# 
# markerGenes= markerGenes[!grepl(pattern='Microglia_',names(markerGenes))]
# 
# markerGenes = markerGenes %>% lapply(function(x){x[x %>% mouse2human %$% humanGene %>% {. %in% rn(trabzuniRegionsExp)}]})
# markerGenes = markerGenes[markerGenes %>% sapply(length) %>% {.>1}]

genes = allCorticalMarkers
genes = genes %>% lapply(function(x){(x %>% mouse2human %$% mouseGene)[x %>% mouse2human %$% humanGene %>% {. %in% rn(trabzuniRegionsExp)}] %>%unique})
genes = genes[genes %>% sapply(length) %>% {.>1}]
genes= genes[!grepl(pattern='Microglia_',names(genes))]
print(genes$`Layer 6b Pyra`)

# white matter grey matter comparison -------------
cortex_white = trabzuniRegionsExp[, c(trabzuniRegionsMeta$brainRegion %in% c('frontal cortex', 'white matter'))]
groups = trabzuniRegionsMeta$brainRegion[trabzuniRegionsMeta$brainRegion %in% c('frontal cortex', 'white matter')]

cortex_white  = data.frame(Gene.Symbol = rownames(cortex_white), cortex_white)

# genes = markerGenes

# names(genes) = translatePublishable(names(genes))
# 
fullEstimate(cortex_white,
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groups,
             outDir='analysis/04.MarkerGeneProfiles/withSingleCells/cortex_white/',
             seekConsensus=FALSE,
             groupRotations=TRUE,
             outlierSampleRemove=FALSE,
             removeNegatives =TRUE,
             comparisons = 'all',
             estimateFile = 'analysis/04.MarkerGeneProfiles/cortex_white//estimations')

#genes$`Pyramidal S100a10` = genes$`Pyramidal S100a10`[genes$`Pyramidal S100a10` %in% allowedProbes]

cortex_whiteEstimate = mgpEstimate(cortex_white,
                                   genes= genes,
                                   geneColName="Gene.Symbol",
                                   groups = groups,
                                   outlierSampleRemove= FALSE,
                                   removeNegatives=TRUE,
                                   seekConsensus = FALSE
)

viable = function(x){
    out = nrow(x)>1
    if(is.null(x)){
        out = FALSE
    }
    return(out)
}

cortex_whiteEstimate$estimates = cortex_whiteEstimate$estimates[(cortex_whiteEstimate$rotations %>% sapply(viable))]
cortex_whiteEstimate$groups = cortex_whiteEstimate$groups[(cortex_whiteEstimate$rotations %>% sapply(viable))]
cortex_whiteEstimate$rotations = cortex_whiteEstimate$rotations[(cortex_whiteEstimate$rotations %>% sapply(viable))]

# white matter gray matter plot -------
# as long as outlier samples are not removed all cortex_whiteEstimate$groups will be the same
frame = data.frame(melt(cortex_whiteEstimate$estimates %>% lapply(scale01)), cortex_whiteEstimate$groups[[1]]) 
names(frame) = c('estimation','cellType','brainRegions')

# get gene counts used in the estimation. rotations list hass rotations of the genes in the PCs. Use it to see
# how many genes were used in the estimation
genesUsed = cortex_whiteEstimate$rotations %>% sapply(nrow)

numberedNames = paste0(names(genesUsed),'\n(n genes = ', genesUsed, ')')
names(numberedNames) = names(genesUsed)
# ogbox toColor replaces elements of a 
frame$cellType = ogbox::replaceElement(frame$cellType,
                                       numberedNames)$newVector

#order = translatePublishable(cellOrder)
# order = cellOrder
#frame %<>% mutate(cellType = cellType %>% 
#                      factor(levels =unique(frame$cellType)[unique(frame$cellType) %>% str_extract('^.*?(?=\n)') %>% match(order,.)] %>% trimNAs ))

# estims = read.table('analysis/04.MarkerGeneProfiles/cortex_white/estimations', header=T,sep='\t')
# estims[1:(ncol(estims)-1)] = apply(estims[1:(ncol(estims)-1)],2,scale01)
# 
# list.

#frame =melt(estims)
#names(frame) = c('brainRegions','cellType','estimation')


frame = frame %>% group_by(cellType)
maxY = frame %>% summarise(max(estimation))
ps= by(frame, frame$cellType, function(x){
    a1 = x %>% filter(brainRegions == 'frontal cortex') %>% ungroup %>% select(estimation) %>% unlist
    a2 = x %>% filter(brainRegions == 'white matter') %>% ungroup %>%  select(estimation) %>% unlist
    wilcox.test(a1,a2)$p.value
})

ps = p.adjust(ps,method='fdr')
markers = ogbox::signifMarker(ps)
signifFrame = data.frame(markers, x= 1.5,y= 1,cellType = names(ps))
# signifFrame %<>% filter(cellType !='PyramidalCorticoThalam')
# frame %<>% 
#     #filter(cellType!='PyramidalCorticoThalam') %>% 
#     ungroup %>% 
#     mutate(cellType = cellType %>% droplevels) %>% group_by(cellType)



p  = frame %>%
    ggplot(aes(x=brainRegions, y = estimation)) + facet_wrap(~cellType,ncol=5) + 
    theme_cowplot(17) + 
    geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    geom_boxplot(width=0.1,fill = 'lightblue') + 
    # geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 14)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10) + 
    coord_cartesian(ylim = c(-0.03, 1.10))  + xlab('') + ylab('MGP estimation')

dir.create('analysis//04.MarkerGeneProfiles/withSingleCells/publishPlot', showWarnings=FALSE,recursive = TRUE)

ggsave(filename='analysis//04.MarkerGeneProfiles/withSingleCells/publishPlot/cortex_WhiteMatterEstimations.png',p,width=11.31,height=7,units='in')



