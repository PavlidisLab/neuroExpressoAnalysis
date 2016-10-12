#' @export
s100a10exception = function(restDir=NULL, genelist = NULL, cores = 1){
    if (detectCores()<cores){ 
        cores = detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }
    registerDoMC(cores)
    
    reactAstro = allowedProbesS100a10
    if (!is.null(restDir)){
        fileNames = list.files(restDir, recursive =T )
        fileNames = fileNames[grepl('S100a10',fileNames)]
        #for(i in fileNames){
        foreach (i = fileNames) %dopar% {
            markerGenes = tryCatch({read.table(paste0(restDir,'/',i))},
                                   error = function(e){
                                       NULL
                                   })
            if(is.null(markerGenes)){
                return()
            }
            markerGenesLeft = markerGenes[markerGenes$V1 %in% reactAstro,]
            write.table(markerGenesLeft, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
        }
    }
    
    # just apply to a single microglia list
    if (!is.null(genelist)){
        return(geneList[geneList %in% reactAstro])
    }
}