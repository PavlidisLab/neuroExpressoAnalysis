#devtools::load_all()
library(homologene)
library(data.table)

markerGenes = mouseMarkerGenes$Cortex
markerGenes= markerGenes[!grepl(pattern='Microglia_',names(markerGenes))]

rnaCoexist(TasicMouseExp,
           tresholds= NULL,
           genes = markerGenes,
           dupResolve = TRUE,
           dataOut = paste0('analysis//02.SingleCellCoexpression/Tasic/Output/'),
           plotOut = paste0('analysis//02.SingleCellCoexpression/Tasic/Plots/'),
           cores  = 16)


rnaCoexist(ZeiselMouseExp,
           tresholds= NULL,
           genes = markerGenes,
           dupResolve = TRUE,
           dataOut = paste0('analysis//02.SingleCellCoexpression/Zeisel/Output/'),
           plotOut = paste0('analysis//02.SingleCellCoexpression/Zeisel/Plots/'),
           cores  = 16)

rnaCoexist(DarmanisHumanExp,
           tresholds= NULL,
           genes = markerGenes %>% lapply(function(x){mouse2human(x) %$% humanGene}),
           dupResolve = TRUE,
           dataOut = paste0('analysis//02.SingleCellCoexpression/Darmanis/Output/'),
           plotOut = paste0('analysis//02.SingleCellCoexpression/Darmanis/Plots/'),
           cores  = 16)

# table generation ------------
zeiselSinglePs = fread('analysis//02.SingleCellCoexpression/Zeisel/Output/realProbs',data.table=F)
tasicSinglePs = fread('analysis//02.SingleCellCoexpression/Tasic//Output/realProbs',data.table=F)
darmanisSinglePs = fread('analysis//02.SingleCellCoexpression/Darmanis/Output/realProbs',data.table=F)

zeiselSinglePs = zeiselSinglePs[-3]
tasicSinglePs = tasicSinglePs[-3]
darmanisSinglePs = darmanisSinglePs[-3]

frame = merge(merge(zeiselSinglePs,tasicSinglePs,by='V1',all=TRUE),darmanisSinglePs,by='V1',all=TRUE)
rownames(frame) = frame$V1
frame = frame[-1]

names(frame) = c('p value','Gene Count', 'p value','Gene Count','p value','Gene Count')
frame %<>% fixTable


frame = data.frame('Cell type' =zeiselSinglePs[1] ,frame,check.names=F)
frame = frame[ogbox::trimNAs(match(cellOrder,frame$V1)),]
frame$V1 = publishableNameDictionary$ShinyNames[match(frame$V1,
                                                      publishableNameDictionary$PyramidalDeep)]

frame %>% write.design('analysis//02.SingleCellCoexpression/singleCellTable.tsv')
