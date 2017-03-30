#' @export
microglialException = function(restDir=NULL, genelist = NULL, cores = 1){
    if (detectCores()<cores){ 
        cores = detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }
    registerDoMC(cores)
    # genes effected by old age and LPS stimulation
    
    effectedGenes = read.table('data-raw/GOAD/LPS-STIMULATED MICROGLIA (4HR) vs. CONTROL MICROGLIA.tsv',
                               header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>%
        filter(Adjusted.p.value<=0.05) %$% Gene.symbol
    
    activationGenes = read.table('data-raw/GOAD/LPS-STIMULATED MICROGLIA (4HR) vs. CONTROL MICROGLIA.tsv',
                                 header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>%
        filter(Adjusted.p.value<=0.05, Fold.Change >1 ) %$% Gene.symbol
    
    deactivationGenes = read.table('data-raw/GOAD/LPS-STIMULATED MICROGLIA (4HR) vs. CONTROL MICROGLIA.tsv',
                                   header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>%
        filter(Adjusted.p.value<=0.05, Fold.Change <1 ) %$% Gene.symbol
    
    
    if (!is.null(restDir)){
        fileNames = list.files(restDir, recursive =T )
        fileNames = fileNames[grepl('Microglia$',fileNames)]
        #for(i in fileNames){
        foreach (i = fileNames) %dopar% {
            micro = read.table(paste0(restDir,'/',i))
            microAll = micro[!toupper(micro$V1) %in% effectedGenes,]
            write.table(microAll, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
            
            actiMicro = micro[toupper(micro$V1) %in% activationGenes,]
            write.table(actiMicro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'_activation'))
            
            deactiMicro = micro[toupper(micro$V1) %in% deactivationGenes,]
            write.table(deactiMicro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'_deactivation'))
            
        }
    }
    
    # just apply to a single microglia list
    if (!is.null(genelist)){
        return(genelist[!toupper(genelist) %in% effectedGenes])
    }
}
