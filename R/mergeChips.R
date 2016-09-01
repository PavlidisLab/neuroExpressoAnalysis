#' Merge two affymetrix chips
#' @description Takes two affy \link[Biobase]{AffyBatch} objects and converts them 
#' into an \link[Biobase]{ExpressionSet} object. It subsets the probes with the same
#' identifiers and discards any probes that are not shared between the two chips.
#' Intended use is merging data from chips and their later versions. 
#' @param affy1 An \link[Biobase]{AffyBatch} object
#' @param affy2 Another \link[Biobase]{AffyBatch} object, sharing probes with aff1
#' @param allowIntersect If TRUE, it will not return an error when one of the chips 
#' do not include all probesets from the other.
#' @export
mergeChips = function(affy1,affy2,allowIntersect=FALSE){
    PNList1 = affy::probeNames(affy1)
    PNList2 = affy::probeNames(affy2)
    
    if (length(PNList2) >length(PNList1)){
        affyTemp = affy2
        affy2 = affy1
        affy1 = affyTemp
        PNList1 = affy::probeNames(affy1)
        PNList2 = affy::probeNames(affy2)
    }
    
    
    subsetList = PNList1[PNList1 %in% PNList2]
    
    
    
    subsetPm = affy::pm(affy1, unique(subsetList))
    subsetPmOldOrdered = affy::pm(affy2, unique(subsetList))
    
    if (!allowIntersect){
        subsetPmOld = affy::pm(affy2)
    } else {
        subsetPmOld = affy::pm(affy2,unique(subsetList))
    }
    
    allPm = cbind(subsetPm, subsetPmOldOrdered)
    
    rownames(allPm) = rownames(subsetPmOld)
    
    subset = NULL
    verbose = TRUE
    destructive = TRUE
    normalize = TRUE
    background = TRUE
    bgversion = 2
    ngenes = length(geneNames(affy2))
    allpNList = split(0:(length(PNList2) - 1), PNList2)
    
    
    #rownames(allPm) = 1:nrow(allPm)
    exprs <- .Call("rma_c_complete", allPm, 
                   allpNList, ngenes, normalize, background, bgversion, 
                   verbose, PACKAGE = "affy")
    
    phenoD = BiocGenerics::combine(phenoData(affy1), phenoData(affy2))
    annot =  BiocGenerics::annotation(affy2)
    protocolD = BiocGenerics::combine(protocolData(affy1), protocolData(affy2))
    experimentD = Biobase::experimentData(affy2)
    
    
    newNormalized = new("ExpressionSet", phenoData = phenoD, annotation = annot, 
                        protocolData = protocolD, experimentData = experimentD, 
                        exprs = exprs)
    return(newNormalized)
    
}