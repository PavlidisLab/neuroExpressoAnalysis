
devtools::load_all()
library(purrr)
library(magrittr)

versionOld = 'Markers_1.0'
versionNew = 'Markers_1.1'
referenceGroup = 'CellTypes'


# oldList = readRDS(glue('analysis/01.SelectGenes/{versionOld}/mouseMarkerGenesCombined.rds'))

oldList = pickMarkersAll(glue('analysis/01.SelectGenes/{versionOld}/Markers_Final/{referenceGroup}'))

# newList = readRDS(glue('analysis/01.SelectGenes/{versionNew}/mouseMarkerGenesCombined.rds'))

newList = pickMarkersAll(glue('analysis/01.SelectGenes/{versionNew}/Markers_Final/{referenceGroup}'))
assertthat::assert_that(all(names(oldList) == names(newList)))

logFile = glue('analysis/01.SelectGenes/comparison_{referenceGroup}-{versionOld}-{versionNew}.md')
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

