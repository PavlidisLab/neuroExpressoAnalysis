
devtools::load_all()
library(purrr)
library(magrittr)

versionOld = 'Markers_1.1'
versionNew = 'Markers_1.2'
referenceGroup = 'CellTypes'


# oldList = readRDS(glue('analysis/01.SelectGenes/{versionOld}/mouseMarkerGenesCombined.rds'))

oldList = pickMarkersAll(glue('analysis/01.SelectGenes/{versionOld}/Markers_Final/{referenceGroup}'))

# newList = readRDS(glue('analysis/01.SelectGenes/{versionNew}/mouseMarkerGenesCombined.rds'))

newList = pickMarkersAll(glue('analysis/01.SelectGenes/{versionNew}/Markers_Final/{referenceGroup}'))
assertthat::assert_that(all(names(oldList) == names(newList)))

logFile = glue('analysis/01.SelectGenes/{versionNew}/comparison_{referenceGroup}-{versionOld}-{versionNew}.md')
file.create(logFile,showWarnings = FALSE)
lapply(names(oldList),function(x){
    print(x)
    oldRegion = oldList[[x]]
    newRegion = newList[[x]]
    assertthat::assert_that(all(names(oldRegion) == names(newRegion)))
    names(oldRegion) %>% lapply(function(y){
        oldGenes = oldRegion[[y]]
        newGenes = newRegion[[y]]
        intersection = intersect(oldGenes,newGenes)
        exclusiveOld = oldGenes[!oldGenes %in% intersection]
        exclusiveNew = newGenes[!newGenes %in% intersection]
        return(list(removed = exclusiveOld,
                    added = exclusiveNew))
    }) -> out
    names(out) = names(oldRegion)
    
    removeCount = out %>% map('removed') %>% sapply(length)
    addCount = out %>% map('added') %>% sapply(length)
    out = out[!(removeCount == 0 & addCount == 0)]
    if(length(out)>0){
        cat('##',x,'\n\n',file = logFile,sep = '',append = TRUE)
        names(out) %>% sapply(function(y){
            cat('###',y,'\n\n',file = logFile,sep = '',append=TRUE)
            cat('Total gene count:',length(newRegion[[y]]),'\n\n',append = TRUE, file = logFile)
            out[[y]] = out[[y]][sapply(out[[y]],length)>0]
            names(out[[y]]) %>% sapply(function(z){
                cat(z,' (',length(out[[y]][[z]]),')',': ',paste(out[[y]][[z]],collapse = ', '),'\n\n',file = logFile,sep = '',append = TRUE)
            })
            
        })
    }
    
    return(out)
}) -> changeList

names(changeList) = names(oldList)

# a shorter repor that gives the number of genes changes

simplerChangeList = changeList %>% unlist(recursive = FALSE) 
simplerNewList= newList %>% unlist(recursive = FALSE)

allCellTypes = names(simplerChangeList) %>% stringr::str_extract('(?<=\\.).*') %>% unique

shortLog = glue('analysis/01.SelectGenes/{versionNew}/comparison_{referenceGroup}-{versionOld}-{versionNew}_short.md')
file.create(shortLog,showWarnings = FALSE)

allCellTypes %>% sapply(function(x){
    cellTypeChange = simplerChangeList[grepl(pattern=x , x= names(simplerChangeList))]
    added = cellTypeChange %>% purrr::map('added') %>% unlist %>% unique
    removed =  cellTypeChange %>% purrr::map('removed') %>% unlist %>% unique
    
    cellTypeGenes = simplerNewList[grepl(pattern=x , x= names(simplerNewList))] %>% unlist %>% unique
    
    # cat('* **',x,'**\n',file = logFile,sep = '',append = TRUE)
    # cat('  - added:',x,'**\n',file = logFile,sep = '',append = TRUE)
    
    return(c('added' = length(added), 'removed' = length(removed),'Total genes in final list' = length(cellTypeGenes)))
}) %>% t %>% knitr::kable(format = 'markdown') %>% cat(file = shortLog,sep='\n')
