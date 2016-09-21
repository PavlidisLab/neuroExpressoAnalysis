#' @export
astrocyteException = function(restDir=NULL, genelist = NULL, cores = 1){
    if (detectCores()<cores){ 
        cores = detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }
    registerDoMC(cores)
    
    frm = data.frame(reactive = c(F,F,F,T,T,T))
    mm = model.matrix(~ reactive,frm)
    fit <- lmFit(astrocytesReactive, mm)
    fit <- eBayes(fit)
    reactAstro = topTable(fit, coef=colnames(fit$design)[2],
                          #lfc = log(1,base=2),
                          number = Inf, 
                          p.value = 0.05
    )
    
    reactAstro$Gene.Symbol = rownames(reactAstro)
    reactAstro %<>% filter(logFC>log(5,base=2))
    reactAstro = reactAstro$Gene.Symbol
    
    if (!is.null(restDir)){
        fileNames = list.files(restDir, recursive =T )
        fileNames = fileNames[!grepl('Astrocyte$',fileNames)]
        #for(i in fileNames){
        foreach (i = fileNames) %dopar% {
            markerGenes = read.table(paste0(restDir,'/',i))
            markerGenesLeft = markerGenes[!toupper(markerGenes$V1) %in% reactAstro,]
            write.table(markerGenesLeft, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
        }
    }
    
    # just apply to a single microglia list
    if (!is.null(genelist)){
        return(geneList[!geneList %in% reactAstro])
    }
}