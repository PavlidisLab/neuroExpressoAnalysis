devtools::load_all()
library(homologene)
library(data.table)
library(ogbox)

# justMicroarray = pickMarkersAll('analysis/01.SelectGenes/FinalGenes1/PyramidalDeep/')
# withSingleCells = pickMarkersAll('analysis/01.SelectGenes/FinalGenes3/PyramidalDeep/')
# withSingleCells = pickMarkersAll('analysis/01.SelectGenes/FinalGenes3/PyramidalDeepNoSingleCellUnlessYouHaveTo//')
# 
# justSingleCells = pickMarkersAll('analysis/01.SelectGenes/QuickJustSingleCell/')

genes = mouseMarkerGenesPyramidalDeep$Cortex
genes = c(genes,list(Pyramidal = mouseMarkerGenes$Cortex$Pyramidal))
#genes = genes %>% lapply(function(x){x[x %>% mouse2human %$% humanGene %>% {. %in% rn(trabzuniRegionsExp)}]})
#genes = genes[genes %>% sapply(length) %>% {.>1}]
genes= genes[!grepl(pattern='Microglia_',names(genes))]

markerGenes=genes

# markerGenes = names(withSingleCells$All) %>% sapply(function(x){
#     c(withSingleCells$All[[x]],justSingleCells$PyramidalDeep[[x]]) %>% unique
# },simplify=FALSE)

markerGenes= markerGenes[!grepl(pattern='Microglia_',names(markerGenes))]
markerGenes = markerGenes[markerGenes %>% sapply(length) %>% {.>1}]

rnaCoexist(TasicMouseExp,
           tresholds= NULL,
           genes = markerGenes,
           dupResolve = TRUE,
           dataOut = paste0('analysis//02.SingleCellCoexpression/WithSingleCells/Tasic/Output/'),
           plotOut = paste0('analysis//02.SingleCellCoexpression/WithSingleCells/Tasic/Plots/'),
           cores  = 15)


rnaCoexist(ZeiselMouseExp,
           tresholds= NULL,
           genes = markerGenes,
           dupResolve = TRUE,
           dataOut = paste0('analysis//02.SingleCellCoexpression/WithSingleCells/Zeisel/Output/'),
           plotOut = paste0('analysis//02.SingleCellCoexpression/WithSingleCells/Zeisel/Plots/'),
           cores  = 15)

markerGenesH = markerGenes %>% lapply(function(x){mouse2human(x) %$% humanGene})
markerGenesH %<>% sapply(function(x){x[x %in% rn(DarmanisHumanExp)]})
markerGenesH = markerGenesH[markerGenesH %>% sapply(length) %>% {.>1}]

rnaCoexist(DarmanisHumanExp,
           tresholds= NULL,
           genes = markerGenesH,
           dupResolve = TRUE,
           dataOut = paste0('analysis//02.SingleCellCoexpression/WithSingleCells/Darmanis/Output/'),
           plotOut = paste0('analysis//02.SingleCellCoexpression/WithSingleCells/Darmanis/Plots/'),
           cores  = 16)

# table generation ------------
zeiselSinglePs = fread('analysis//02.SingleCellCoexpression/WithSingleCells/Zeisel/Output/realProbs',data.table=F)
tasicSinglePs = fread('analysis//02.SingleCellCoexpression/WithSingleCells/Tasic//Output/realProbs',data.table=F)
darmanisSinglePs = fread('analysis//02.SingleCellCoexpression/WithSingleCells/Darmanis/Output/realProbs',data.table=F)

zeiselSinglePs = zeiselSinglePs[-3]
tasicSinglePs = tasicSinglePs[-3]
darmanisSinglePs = darmanisSinglePs[-3]

frame = merge(merge(zeiselSinglePs,tasicSinglePs,by='V1',all=TRUE),darmanisSinglePs,by='V1',all=TRUE)
rownames(frame) = frame$V1
frame = frame[-1]

names(frame) = c('p value','Gene Count', 'p value','Gene Count','p value','Gene Count')
frame %<>% fixTable


frame = data.frame('Cell type' =zeiselSinglePs[1] ,frame,check.names=F)
#frame = frame[ogbox::trimNAs(match(cellOrder,frame$V1)),]
# frame$V1 = publishableNameDictionary$ShinyNames[match(frame$V1,
#                                                       publishableNameDictionary$PyramidalDeep)]

frame %>% write.design('analysis//02.SingleCellCoexpression/WithSingleCells/singleCellTable.tsv')

