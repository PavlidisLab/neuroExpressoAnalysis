# this version merges cell types into a single sample then looks for a variable.
# specific to our cell type data. not ideal but oh well...
#' @export
mostVariableCT = function(whichFile,outFile=NULL,cellTypeColumn, design){
    if (is.character(whichFile)){
        allDataPre = ogbox::read.exp(whichFile)
    } else{
        allDataPre = whichFile
    }
    
    if (is.character(design)){
        design = ogbox::read.design(design)
    } else{
        design = design
    }

    
    list[,exprData]= sepExpr(allDataPre)
    
    cellTypes = trimNAs(unique(design[,cellTypeColumn]))
    
    cellTypeExpr = lapply(cellTypes,function(x){
        apply(exprData[,design[,cellTypeColumn] %in% x,drop=F],1,mean)
    })
    exprData = as.data.frame(cellTypeExpr)
    # remove the ones with highest expression below 6
    rowmax = apply(exprData, 1, max)
    discludeGenes = (rowmax<6)
    allDataPre = allDataPre[!discludeGenes,]
    exprData = exprData[!discludeGenes,]
    
    # ignore multiple matching probesets while mostVariable selection
    allDataMulti = allDataPre[grepl('[|]',allDataPre$Gene.Symbol),]
    exprData = exprData[!grepl('[|]',allDataPre$Gene.Symbol),]
    allDataPre = allDataPre[!grepl('[|]',allDataPre$Gene.Symbol),]
    
    # you bloody idiot... taken from lila
    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    allDataPre = allDataPre[decreasingVar,]
    allDataPre = allDataPre[!duplicated(allDataPre$Gene.Symbol),]
    
    # add the multiple matching probesets back
    allDataPre = rbind(allDataPre,allDataMulti)
    allDataPre = allDataPre[!allDataPre$Gene.Symbol=='',]
    if(!is.null(outFile)){
        write.csv(allDataPre, file = outFile, row.names=FALSE)
    }
    invisible(allDataPre)
}

# this function is a generic function that looks for the most variable probeset
# of a gene. unlike the previous one, it takes in objects and outputs objects 
#' @export
mostVariable = function(allDataPre,genes = 'Gene.Symbol', threshold = 6, threshFun = max){
    list[,exprData]= sepExpr(allDataPre)
    rowmax = apply(exprData, 1, threshFun)
    discludeGenes = (rowmax<threshold)
    allDataPre = allDataPre[!discludeGenes,]
    exprData = exprData[!discludeGenes,]
    
    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    allDataPre = allDataPre[decreasingVar,]
    if (class(allDataPre)[1]=='data.table'){
        allDataPre = allDataPre[!duplicated(allDataPre[,genes, with=F]),]
    } else {
        allDataPre = allDataPre[!duplicated(allDataPre[,genes]),]
        
    }
    allDataPre = allDataPre[!allDataPre[,genes]=='',]
    return(allDataPre)
}