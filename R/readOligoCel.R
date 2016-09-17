#' @export
readOligoCel = function(GSMs, gemmAnnot,fileOut=NULL, celdir){
    cels = oligoClasses::list.celfiles(celdir,listGzipped=T)
    whichCels = sapply(GSMs, function(x){which(grepl(x,cels))})
    sampleCels = cels[whichCels]
    affyRaw = oligo::read.celfiles(paste0(celdir,'/',sampleCels))
    exonTS <- oligo::rma(affyRaw, target = "core")
    rm(affyRaw)
    featureData(exonTS) <- getNetAffx(exonTS, "transcript")
    #View the features of the obtained data
    # exonTS
    #Extract the expression data
    
    aned = gemmaAnnotOligo(normalized=exonTS,chipFile=gemmaAnnot)
    
    #Create the expression file
    if (!is.null(fileOut)){
        write.csv(aned, fileOut, row.names=FALSE)
    }
    invisible(aned)
}