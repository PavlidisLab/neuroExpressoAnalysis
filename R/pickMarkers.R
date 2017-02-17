#' @export
pickMarkersAll = function(genesLoc,lilah=F,regex='*'){
    allGenLocs = list.dirs(genesLoc)
    allGenLocs = allGenLocs[-1]
    allGenLocs = grep(regex,allGenLocs,value=T)
    geneLists = lapply(allGenLocs, pickMarkers)
    names(geneLists) = basename(allGenLocs)
    return(geneLists)
}

#' Pick marker genes out of candidate lists
#'
#' Picks marker genes out of candidate lists. Behaves differently based on the number of columns the files have.
#' @export
pickMarkers = function(geneLoc, rotationThresh = 0.95,silhouette = 0.5,foldChange = 10,lilah = F){
    filenames = list.files(geneLoc,include.dirs = FALSE)
    filenames = filenames[!filenames %in% 'removed']
    fileContents = lapply(paste0(geneLoc,'/', filenames), function(x){
        tryCatch(
            read.table(x,stringsAsFactors=FALSE),
            error = function(e){
                NULL
            })
    })
    names(fileContents) = filenames

    # empty first file protection
    lengths = sapply(fileContents,len)
    destinedToBeFirst = which.max(lengths>0)

    theFirst = fileContents[1]
    fileContents[1] = fileContents[destinedToBeFirst]
    fileContents[destinedToBeFirst] = theFirst

    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = names(fileContents)

    if (ncol(fileContents[[1]])==3 & lilah == F){
        # this if for a combination of fold change and silhouette coefficient
        for (i in 1:length(fileContents)){
            tempContent = fileContents[[i]][!fileContents[[i]][,1] %in%
                                                unlist(sapply((1:len(fileContents))[-i], function(x){
                                                    fileContents[[x]][,1]
                                                })),]

            geneList[[i]] = as.character(tempContent$V1[as.numeric(as.character(tempContent$V3))>silhouette
                                                        & as.numeric(as.character(tempContent$V2))>log(foldChange,base=2)])
        }
    }else if (ncol(fileContents[[1]])==3 & lilah == T){
        # this if for lilah's selection method
        for (i in 1:length(fileContents)){
            tempContent = fileContents[[i]][!fileContents[[i]][,1] %in%
                                                unlist(sapply((1:len(fileContents))[-i], function(x){
                                                    fileContents[[x]][,1]
                                                })),]
            geneList[[i]] = as.character(tempContent$V1[as.numeric(as.character(tempContent$V3))*
                                                            as.numeric(as.character(tempContent$V2))>2])
        }
    } else if (ncol(fileContents[[1]])==1){
        # this is for a mere gene list
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else if(ncol(fileContents[[1]])==2){
        # this is for selection of percentages from confidence output
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>rotationThresh])
        }
    } else if(ncol(fileContents[[1]])==4){
        # this is for 10 fold changed markers. none of the other collumns matter. just get the genes dammit
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else {
        stop('What kind of gibberish is this')
    }

    puristList = vector(mode = 'list', length = length(geneList))

    # if the file only has a single column, this means it is the final list. consider putting this to a new
    # place later since the point of this function was supposed to be filtering the genes as well
    if (ncol(fileContents[[1]])!=1){
        for (i in 1:length(geneList)){
            puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
        }
    } else {
        puristList = geneList
    }
    names(puristList) = names(geneList)
    puristList = lapply(puristList, as.character)

    theFirst = fileContents[1]
    fileContents[1] = fileContents[destinedToBeFirst]
    fileContents[destinedToBeFirst] = theFirst


    return(puristList)
}

# lapply(c('Eno2','Mog'), findGene, list)
#' Find if a marg
#' @export
findGene = function(gene,list){
    out = lapply(list, function(x){
        findInList(gene,x)
    })
    matches = out[lapply(out,len)>0]
    if (len(matches)<1){
        return(NULL)
    }
    matches = sapply(1:len(matches), function(i){
        paste0(names(matches[i]),'_', names(list[[names(matches[i])]][matches[[i]]]))
    })
    return(matches)
}

