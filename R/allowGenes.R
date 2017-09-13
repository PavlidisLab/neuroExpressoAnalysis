#'@export
allowGenes = function(restDir=NULL, genelist = NULL, allowedGenes, regex = NULL,  cores = 1){
    if(!is.na(detectCores())){
        if (detectCores()<cores){ 
            cores = detectCores()
            print('max cores exceeded')
            print(paste('set core no to',cores))
        }
    }
    cl<-parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    
    if (!is.null(restDir)){
        fileNames = list.files(restDir, recursive =T )
        if(!is.null(regex)){
            fileNames = fileNames[grepl(regex,fileNames)]
        }
        #for(i in fileNames){
        foreach (i = fileNames) %dopar% {
            print(i)
            markerGenes = tryCatch({read.table(paste0(restDir,'/',i),stringsAsFactors=FALSE)},
                                   error = function(e){
                                       NULL
                                   })
            if(is.null(markerGenes)){
                return()
            }
            markerGenesLeft = markerGenes[markerGenes$V1 %in% allowedGenes,]
            write.table(markerGenesLeft, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
        }
    }
    
    # just apply to a single microglia list
    if (!is.null(genelist)){
        return(geneList[geneList %in% allowedGenes])
    }
}
