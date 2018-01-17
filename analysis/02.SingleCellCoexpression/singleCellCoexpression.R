devtools::load_all()
library(homologene)
library(data.table)

markerGenes = mouseMarkerGenesCombined$Cortex
markerGenes= markerGenes[!grepl(pattern='Microglia_',names(markerGenes))]

markerGenes = markerGenes[!grepl(pattern = '(?!^Pyramidal$)Pyra',x = names(markerGenes),perl = TRUE)]


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


frame = data.frame('Cell type' =rownames(frame) ,frame,check.names=F)
frame = frame[ogbox::trimNAs(match(cellOrder,frame$`Cell type`)),]
frame$`Cell type` = publishableNameDictionary$ShinyNames[match(frame$`Cell type`,
                                                      publishableNameDictionary$PyramidalDeep)]

frame %>% write.design('analysis//02.SingleCellCoexpression/singleCellTable.tsv')
