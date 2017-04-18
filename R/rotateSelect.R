
#' @export
rotateSelect = function(rotationOut,rotSelOut,cores=4, lilah=F, ...){
    
    # so that I wont fry my laptop
    if(!is.na(detectCores())){
        if (detectCores()<cores){ 
            cores = detectCores()
            print('max cores exceeded')
            print(paste('set core no to',cores))
        }
    }
    registerDoMC(cores)
    
    dirFols = list.dirs(rotationOut, recursive = F)
    
    loopAround = list.dirs(dirFols[1],full.names = F)
    
    # in server version full.names input of list.dirs do not work. This fixes it. Might add this to ogbox as an overwrite.
    loopAround = gsub(paste0(dirFols[1],'/'),'',loopAround)
    loopAround = loopAround[!grepl(dirFols[1],loopAround)]
    
    loopAround = loopAround[!loopAround %in% '']
    #loopAround = loopAround [-which(loopAround %in% c('Relax','Marker',''))]
    #loopAroundRel = loopAround[grepl('Relax',loopAround)]
    #loopAroundMar = loopAround[grepl('Marker',loopAround)]
    dir.create(paste0(rotSelOut), showWarnings = F)
    
    # for relaxed selection. forces unique selection with matching criteria in puristOut
    # for (i in loopAroundRel){
    foreach (i  = loopAround) %dopar% {
        print(i)
        dir.create(paste0(rotSelOut,'/',i),recursive = T, showWarnings = F)
        files = list.files(paste0(dirFols[1],'/',i))
        # remove the list of removed samples from the mix
        files = files[!files %in% 'removed']
        
        pureConfidence = vector(mode = 'list', length =len(files))
        for (j in dirFols){
            #print(paste0(j,'/',i))
            pureSample = pickMarkers(paste0(j,'/',i), lilah, 
                                   ...
            )
            pureConfidence = mapply(c,pureSample,pureConfidence,SIMPLIFY=FALSE)
            #print(names(pureConfidence))
            if(len(pureConfidence)>len(files)){
                stop('dafaq man')
            }
        }
        confidence = lapply(pureConfidence,function(x){table(x)/len(dirFols)})
        
        #         if (any(grepl('removed',names(confidence)))){
        #             confidence = confidence[-which(grepl('removed',names(confidence)))]
        #         }
        
        for (j in 1:len(confidence)){
            # genes = names(confidence[[j]])[confidence[[j]]>0.95]
            write.table(confidence[[j]],row.names=F,quote=F,col.names=F,
                        file = paste0(rotSelOut,'/',i,'/',
                                      names(confidence)[j]))
        }
        #return(invisible())
    }
    
    # for marker genes. just looks at the list and frequencies
    # for (i in loopAroundMar){
#     foreach (i = loopAroundMar) %dopar%{
#         print(i)
#         
#         dir.create(paste0(rotSelOut,'/',i),recursive = T, showWarnings = F)
#         files = list.files(paste0(dirFols[1],'/',i))
#         fileOuts = vector(mode = 'list', length = len(files))
#         for (j in files){
#             genes = vector()
#             for (k in dirFols){
#                 daFile = tryCatch({read.table(paste0(k,'/',i,'/',j) ,header=F)},
#                                   error=function(e){
#                                       daFile=data.frame()
#                                   })
#                 genes = c(genes, as.character(daFile$V1))
#             }
#             
#             geneCounts = table(genes)
#             confidence = geneCounts/len(dirFols)
#             
#             write.table(as.df(confidence), file = paste0(rotSelOut,'/',i,'/',j),
#                         row.names = F, col.names = F, quote=F, sep='\t')
#         }
#         return(invisible())
#         
#     }
    
    
}