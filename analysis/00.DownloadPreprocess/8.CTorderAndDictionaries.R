
firstChip = TRUE
singleCell = TRUE
typeSets = c('PyramidalDeep','CellTypes')
# now load the data and place it in the package
n_expressoSamples =  ogbox::read.design('data-raw/Mouse_Cell_Type_Data/n_expressoSamples.tsv')

TasicSamples = ogbox::read.design('data-raw/Mouse_Cell_Type_Data/meltedSingleCells.tsv')

samples = rbind(n_expressoSamples,TasicSamples)

# order of cell types that makes sense
cellOrder = samples %>%
    arrange(MajorType,Neurotransmitter,PyramidalDeep) %>% 
    filter(!is.na(PyramidalDeep)) %>% .$PyramidalDeep %>% unique
cellOrder = c(cellOrder[1:27],'Pyramidal',cellOrder[28:36])
#cellOrder = c(cellOrder[1],'AstrocyteReactive','AstrocyteInactive', cellOrder[2:len(cellOrder)])
devtools::use_data(cellOrder,overwrite=TRUE)

publishableNameDictionary = samples %>% filter(!is.na(PyramidalDeep)) %>%  select(PyramidalDeep, ShinyNames) %>% unique
publishableNameDictionary = rbind(publishableNameDictionary,c('Pyramidal','Pyramidal'))

devtools::use_data(publishableNameDictionary,overwrite=TRUE)


# select simple markers for verification
if(firstChip & singleCell){
    for (x in typeSets){
        # for neuroExpresso
        cortex = memoReg(n_expressoSamples,regionNames = 'Region',groupNames = x,regionHierarchy = regionHierarchy)$Cortex
        
        n_Exp = n_expressoExpr %>% filter(!grepl('\\|',Gene.Symbol))
        list[gene, exp] = sepExpr(n_Exp)
        rownames(exp) = gene$Gene.Symbol
        exp = exp[!is.na(cortex)]
        n_samples = n_expressoSamples[!is.na(cortex),]
        NeuroExpressoPrimaryMean = n_samples[[x]] %>% unique %>% trimNAs %>% lapply(function(y){
            exp[, n_samples[[x]] %in% y] %>% apply(1,mean)
        }) %>% as.data.frame
        names(NeuroExpressoPrimaryMean) =  n_samples[[x]] %>% unique %>% trimNAs 
        rownames(NeuroExpressoPrimaryMean) = gene$Gene.Symbol
        # use_data(NeuroExpressoPrimaryMean,overwrite = TRUE)
        
        nxSimpleMarkers = NeuroExpressoPrimaryMean %>% apply(1,which.max) %>% names(NeuroExpressoPrimaryMean)[.]
        names(nxSimpleMarkers) = rn(NeuroExpressoPrimaryMean)
        teval(paste0('nxSimpleMarkers_',x,'<<-nxSimpleMarkers'))
        teval(paste0('use_data(nxSimpleMarkers_',x,",overwrite=TRUE)"))
        
        # for Tasic
        singleCells = ogbox::read.design('data-raw/Mouse_Cell_Type_Data/singleCellMatchings.tsv')
        
        tasicCellTypeMeans = singleCells[[x]] %>% unique %>% lapply(function(y){
            cluster = (singleCells$Tasic[singleCells[[x]] %in% y]) %>% str_split(', ') %>% {.[[1]]}
            TasicPrimaryMean[cluster] %>% apply(1,mean) 
        }) %>% as.data.frame()
        names(tasicCellTypeMeans) =  singleCells[[x]] %>% unique %>% trimNAs()
        tasicSimpleMarkers = tasicCellTypeMeans %>% apply(1,which.max) %>% names(tasicCellTypeMeans)[.]
        names(tasicSimpleMarkers) = rn(tasicCellTypeMeans)
        teval(paste0('tasicSimpleMarkers_',x,'<<-tasicSimpleMarkers'))
        teval(paste0('use_data(tasicSimpleMarkers_',x,",overwrite=TRUE)"))
        
    }
}
