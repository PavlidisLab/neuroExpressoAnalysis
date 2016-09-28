#' @export
banGenes = function(restDir=NULL, genelist = NULL, bannedGenes, cores = 1){
    if (detectCores()<cores){ 
        cores = detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }
    registerDoMC(cores)
        
    if (!is.null(restDir)){
        fileNames = list.files(restDir, recursive =T )
        #for(i in fileNames){
        foreach (i = fileNames) %dopar% {
            markerGenes = tryCatch({read.table(paste0(restDir,'/',i))},
                                   error = function(e){
                                       NULL
                                   })
            if(is.null(markerGenes)){
                return()
            }
            markerGenes = read.table(paste0(restDir,'/',i))
            markerGenesLeft = markerGenes[!toupper(markerGenes$V1) %in% bannedGenes,]
            write.table(markerGenesLeft, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
        }
    }
    
    # just apply to a single microglia list
    if (!is.null(genelist)){
        return(geneList[!geneList %in% bannedGenes])
    }
}