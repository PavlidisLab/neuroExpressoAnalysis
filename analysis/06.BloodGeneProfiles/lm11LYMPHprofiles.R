library(ogbox)
devtools::load_all()
library(homologene)




lm11GenesHuman = humanBloodMarkers$lm11
lm11GenesHuman = lm11GenesHuman[c("B Cells",
                                  'CD4',
                                  'CD8 T cells')]


names(lm11GenesHuman)[2] = 'CD4 T cells'
names(lm11GenesHuman)[1] = 'B cells'

lm11GenesMouse = mouseBloodMarkers$lm11
lm11GenesMouse = lm11GenesMouse[c("B Cells",
                                  'CD4',
                                  'CD8 T cells')]
names(lm11GenesMouse)[2] = 'CD4 T cells'
names(lm11GenesMouse)[1] = 'B cells'


geneLists = list(lm11GenesHuman,
                 lm11GenesMouse %>% lapply(function(x){mouse2human(x)%$% humanGene}))

names(geneLists) =  c('lm11GenesHumanLYMPH',
                      'lm11GenesMouseLYMPH')
# renaming for matching purposes


LYMPHCibersort %<>% arrange(Cell.type, Sample.ID)

cors = lapply(1:len(geneLists), function(i){
    dir.create(paste0('analysis//06.BloodGeneProfiles/plots/',names(geneLists)[i]),recursive=TRUE,showWarnings=FALSE)
    superImpose(geneLists[[i]],
                expression = LYMPHexpr,
                CIBERSORTpred = LYMPHCibersort[],
                geneSymbol = 'Gene.Symbol',
                outDir = paste0('analysis/06.BloodGeneProfiles/plots/',names(geneLists)[i],'/'))
})

dir.create('analysis//06.BloodGeneProfiles/publishTable',showWarnings=FALSE)
cors[[2]] = cors[[2]][,-1,drop=FALSE]
write.table(cors, 'analysis//06.BloodGeneProfiles/publishTable/LYMPHLM11.tsv',sep='\t', quote= FALSE,row.names=TRUE)
