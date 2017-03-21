# now load the data and place it in the package
n_expressoExpr = ogbox::read.exp('data-raw//Mouse_Cell_Type_Data//n_expressoExpr.csv')
n_expressoSamples =  ogbox::read.design('data-raw/Mouse_Cell_Type_Data/n_expressoSamples.tsv')

TasicSamples = ogbox::read.design('data-raw/Mouse_Cell_Type_Data/meltedSingleCells.tsv')

samples = rbind(n_expressoSamples,TasicSamples)

# order of cell types that makes sense
cellOrder = samples %>%
    arrange(MajorType,Neurotransmitter,PyramidalDeep) %>% 
    filter(!is.na(PyramidalDeep)) %>% .$PyramidalDeep %>% unique
#cellOrder = c(cellOrder[1],'AstrocyteReactive','AstrocyteInactive', cellOrder[2:len(cellOrder)])
devtools::use_data(cellOrder,overwrite=TRUE)

publishableNameDictionary = samples %>% filter(!is.na(PyramidalDeep)) %>%  select(PyramidalDeep, ShinyNames) %>% unique
publishableNameDictionary = rbind(publishableNameDictionary,c('Pyramidal','Pyramidal'))

devtools::use_data(publishableNameDictionary,overwrite=TRUE)