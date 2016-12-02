#' @export
bioGender <- function(aned,geneColName = 'Gene.Symbol', probeColName = 'Probe', maleGenes = c('RPS4Y1','KDM5D') , femaleGenes = 'XIST'){
    
    probeM <- aned[[probeColName]][grep(ogbox::regexMerge(maleGenes), aned[[geneColName]])] %>% as.character
    probeF <- aned[[probeColName]][grep(ogbox::regexMerge(femaleGenes), aned[[geneColName]])] %>% as.character
    print(paste('male probes found:', paste(probeM,collapse =', ')))
    print(paste('female probes found:', paste(probeF,collapse =', ')))
    
    
    if(length(probeM) > 0 | length(probeF) > 0){
        dataSex <- t(aned[aned[[probeColName]] %in% c(probeM, probeF),] %>%
                         (ogbox::sepExpr) %>% {.[[2]]})
        
        colnames(dataSex) <- aned[[probeColName]][aned[[probeColName]] %in%  c(probeM, probeF)]
        Clusters <- kmeans(dataSex, centers=2)
        
        Centers <- Clusters$center
        if(length(probeM)>0 & length(probeF)>0){
            if(mean(Centers[1,probeM]) > mean(Centers[1,probeF]) & mean(Centers[2,probeM]) < mean(Centers[2,probeF])){
                Clusters$cluster[Clusters$cluster==1] <- "M"
                Clusters$cluster[Clusters$cluster==2] <- "F"
            } else if(mean(Centers[1,probeM]) < mean(Centers[1,probeF]) & mean(Centers[2,probeM]) > mean(Centers[2,probeF])){
                Clusters$cluster[Clusters$cluster==1] <- "F"
                Clusters$cluster[Clusters$cluster==2] <- "M"
            } else {
                stop("Gender genes disagree, cannot decide about biological gender")
            }
        } else if(length(probeF) == 0) {
            if (mean(Centers[1,]) > mean(Centers[2,])) {
                Clusters$cluster[Clusters$cluster==1] <- "M"
                Clusters$cluster[Clusters$cluster==2] <- "F"
            } else {
                Clusters$cluster[Clusters$cluster==1] <- "F"
                Clusters$cluster[Clusters$cluster==2] <- "M"
            }  
        } else if (length(probeM) == 0) {
            if (mean(Centers[1,]) > mean(Centers[2,])){
                Clusters$cluster[Clusters$cluster==1] <- "F"
                Clusters$cluster[Clusters$cluster==2] <- "M"
            } else {
                Clusters$cluster[Clusters$cluster==1] <- "M"
                Clusters$cluster[Clusters$cluster==2] <- "F"
            }
        } 
        
        BioGender <- Clusters$cluster
        return(BioGender)
    } else {
        print("No sex specific genes on the platform")
    }
}