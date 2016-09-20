library(ogbox)
devtools::load_all()
library(homologene)




lm22GenesHuman = humanBloodMarkers$lm22
lm22GenesHuman = lm22GenesHuman[c('Naïve B cells',
                                'Memory B cells',
                                'CD8 T cells',
                                'CD4 naïve T cells',
                                'CD4 memory T cells-',
                                'CD4 memory T cells+',
                                'NK cells-',
                                'Monos')]
names(lm22GenesHuman)[8] = 'Monocytes'

lm22GenesMouse = mouseBloodMarkers$lm22
lm22GenesMouse = lm22GenesMouse[c('Naïve B cells',
                                 'Memory B cells',
                                 'CD8 T cells',
                                 'CD4 naïve T cells',
                                 'CD4 memory T cells-',
                                 'CD4 memory T cells+',
                                 'NK cells-',
                                 'Monos')]
names(lm22GenesMouse)[8] = 'Monocytes'

geneLists = list(lm22GenesHuman,
                 lm22GenesMouse %>% lapply(function(x){mouse2human(x)%$% humanGene}))
names(geneLists) =  c('lm22GenesHumanPBMC',
                      'lm22GenesMousePBMC')

# renaming for matching purposes
PBMCcibersort$Cell.type[PBMCcibersort$Cell.type == 'Naïve CD4 T cells'] = 'CD4 naïve T cells'
PBMCcibersort$Cell.type[PBMCcibersort$Cell.type == 'Resting memory CD4 T cells'] = 'CD4 memory T cells-'
PBMCcibersort$Cell.type[PBMCcibersort$Cell.type == 'NK cells'] = 'NK cells-'
PBMCcibersort$Cell.type[PBMCcibersort$Cell.type == 'Activated memory CD4 T cells'] = 'CD4 memory T cells+'


cors = lapply(1:len(geneLists), function(i){
    dir.create(paste0('analysis//06.BloodGeneProfiles/plots/',names(geneLists)[i]),recursive=TRUE,showWarnings=FALSE)
    superImpose(geneLists[[i]],
                expression = PBMCexpr,
                CIBERSORTpred = PBMCcibersort,
                geneSymbol = 'GeneSym',
                outDir = paste0('analysis/06.BloodGeneProfiles/plots/',names(geneLists)[i],'/'))
})

dir.create('analysis//06.BloodGeneProfiles/publishTable',showWarnings=FALSE)

write.table(cors, 'analysis//06.BloodGeneProfiles/publishTable/PBMCLM22.tsv',sep='\t', quote= FALSE,row.names=TRUE)
