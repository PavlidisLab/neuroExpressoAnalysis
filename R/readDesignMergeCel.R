# Reads CEL files from a given directory according to a design file
#essential parts of the file are
#1. the first collumn must include the GSMs
#2. platform must be written to the Platform collumn



readMouseCel = function(GSMs,mouseDir='data-raw/cel',gemmaDir = 'data-raw/GemmaAnnots', file=NA){
    allCells = affy::list.celfiles(mouseDir, recursive = T)
    files = sapply(paste0('\\Q',GSMs,'\\E'),function(x){grep(paste0(x,'[.](C|c)(E|e)(L|l)'),allCells,value=T)})
    # notice if files are missing. send out a warning
    if (files %>% sapply(len)  %>% is_greater_than(0) %>% not %>%  any){
        stop('some cel files are missing: ',
             paste((files %>% sapply(len)  %>% is_greater_than(0) %>% not %>% which %>% names),collapse=', '))
    }
    
    
    platforms = unique(dirname(unlist(files)))
    affies =  lapply( rep("AffyBatch", len(platforms)), new )
    
    for (i in 1:len(affies)){
        affies[i] = ReadAffy(filenames = paste0(mouseDir,'/',grep(platforms[i],files, value=T)))
    }
    
    if (len(affies)>1){newNormalized = mergeChips(affies[[1]],affies[[2]])
    } else {newNormalized = affy::rma(affies[[1]])}
    
    # if any one of them are from the old platform, always use the annotation for that one
    if('GPL339' %in% platforms){
        aned = gemmaAnnot(newNormalized, paste0(gemmaDir,'/','GPL339'))
    } else{
        aned = gemmaAnnot(newNormalized, paste0(gemmaDir,'/','GPL1261'))
    }
    aned = aned[!aned$Gene.Symbol == '',]
    
    names(aned) = gsubMult(c('[.](C|c)(E|e)(L|l)','[.](G|g)(z|Z)'),c('',''),names(aned))
    
    if (!is.na(file)){
        write.csv(aned, file, row.names=FALSE)
    }
    invisible(aned)
    
}

#' Read the cell type study metadata file and create the RMA normalized expression matrix.
#' @description This function reads the study design file 'data-raw/Mouse_Cell_Type_Data/cellTypeStudies.tsv' 
#' included in the package and cel files pointed in the file to create the RMA normalized
#'  expression matrix. The products of this function are already loaded with the package. 
#' Use only for replication or upon editing cellTypeStudies.tsv for personal use.
#' @param desFile Address of the study metadata file
#' @param gsm column name that has the GSM identifiers of the samples. The identifiers must
#' be sperated by commas
#' @param normalize The column name that marks which studies should be included.
#' @param celRegex Regular expression that picks the sample names from the gsm column. 
#' Default option is for the cellTypeStudies.tsv. Any modifications to the file should not
#' require changes to the default as long as they are samples from GEO with GSM identifiers.
#' If files that are named differently are added, modify the regex to recognize them.
#' @param expFile name of the output file for gene expression. The file will be a csv
#' @param desOut name of the output file for sample metadata. The file will be a tsv
#' @param gemmaDir The directory that houses the gemma annotations that are created by \link[ogbox]{getGemmaAnnot}
#' @export
readDesignMergeCel = function (desFile = 'data-raw/Mouse_Cell_Type_Data/cellTypeStudies.tsv', 
                               gsm = 'samples',
                               normalize,
                               celRegex = "(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))",
                               celDir ='data-raw/cel/',
                               expFile,
                               desOut, 
                               gemmaDir = 'data-raw/GemmaAnnots'){
    #always have gsms in the first collumn
    
    design = read.table(desFile,quote='',header=T,sep='\t')
    
    design = design[(design[,normalize]),]
    
    gsms = regmatches(design[, gsm], gregexpr(celRegex, design[, gsm],perl=T))
    if (any(lapply(gsms,len)==0)){
        print('No cel names could be captured from the following rows')
        print((1:len(gsms))[lapply(gsms,len)==0])
    }
    
    aned = readMouseCel(unlist(gsms),celDir, gemmaDir)
    
    
    dir.create(dirname(expFile), recursive = TRUE,showWarnings=FALSE)
    dir.create(dirname(desOut), recursive = TRUE,showWarnings=FALSE)
    
    list[genes,exp] = sepExpr(aned)
    #boxplot(aned[,4:ncol(aned)])
    
    header = names(exp)
    
    indexes = vector()
    for (i in 1:length(header)){
        indexes = c(indexes, findInList(header[i], gsms))
    }
    header[!header %in% unlist(gsms)]
    newDesign = data.frame(sampleName = header, originalIndex = indexes, design[indexes,])
    colnames(newDesign) = c('sampleName','originalIndex',names(design))
    #newDesign$originalIndex = as.numeric(newDesign$originalIndex)
    
    #newDesign$age = gsub('~', '', newDesign$age)
    #newDesign$age = gsub('P', '', newDesign$age)
    #newDesign$age = gsub('7-8', '7.5', newDesign$age)
    #newDesign$age[grepl('(precise age not given)',newDesign$age)] = 60
    #newDesign$age = as.numeric(newDesign$age)
    newDesign = newDesign[order(as.numeric(rownames(newDesign))),]
    
    write.table(newDesign, desOut, row.names=FALSE,sep = '\t', quote=F)
    
    #list[genes,expr]=sepExpr(aned)
    
    #expr = expr[match(make.names(design$sampleName)]
    
    exp = exp[,match(make.names(newDesign$sampleName),make.names(colnames(exp)))]
    aned = cbind(genes,exp)
    write.csv(aned, expFile, row.names=FALSE)
    
}


# for changes in original design file that only involves naming. do not run the whole thing again
#'@export
meltDesign = function(desFile, gsm = 'samples', 
                      normalize, celRegex, exprFile, outFile){
    expr = read.exp(exprFile , header = T)
    list[gene,expres] = sepExpr(expr)
    header =  colnames(expres)
    
    design = read.table(desFile,quote='',header=T,sep='\t')
    design = design[(design[,normalize]),]
    gsms = regmatches(design[, gsm], gregexpr(celRegex, design[, gsm],perl=T))
    if (any(lapply(gsms,len)==0)){
        print('No cel names could be captured from the following rows')
        print((1:len(gsms))[lapply(gsms,len)==0])
    }
    
    indexes = vector()
    for (i in 1:length(header)){
        indexes = c(indexes, findInList(header[i], gsms))
    }

    
    if (len(header[!header %in% unlist(gsms)])>0){
        print("I can't find some GSMs in your design file")
        print(header[!header %in% make.names(unlist(gsms))])
    }
    
    newDesign = data.frame(sampleName = header, originalIndex = indexes, design[indexes,])
    
    newDesign = newDesign[match(make.names(colnames(expres)),make.names(newDesign$sampleName)),]
    colnames(newDesign) = c('sampleName','originalIndex',names(design))
    #  newDesign = newDesign[order(as.numeric(rownames(newDesign))),]
    write.table(newDesign, paste0(outFile), row.names=FALSE,sep = '\t', quote=F)
    
}
