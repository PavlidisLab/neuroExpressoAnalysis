load_all()

library(dplyr)
library(viridis)
library(reshape2)
library(ogbox)
library(scales)
library(gplots)
library(stringr)
library(magrittr)
library(pheatmap)
# expression of top 5 heatmap --------

genes = mouseMarkerGenes

genes %<>% lapply(function(x){x[!grepl(pattern='Microglia_',names(x))]})
geneNames = names(genes)
genes = lapply(1:len(genes), function(i){
    nm = names(genes[[i]])
    lapply(1:len(genes[[i]]), function(j){
        print(paste(i,j))
        component1 = names(genes[i])
        component2 = names(genes[[i]][j])
        component1 = paste0(component1,'_')
        if (component1=='All_'){
            component1=''
        }
        scores = read.table(paste0('analysis/01.SelectGenes/Quick/', component1,'PyramidalDeep/',component2))
        (scores %>% filter(V1 %in% genes[[i]][[j]]))[1:5,]$V1 %>% unlist %>% as.char
    }) ->out
    names(out) = nm
    return(out)
})

names(genes) = geneNames
genes = lapply(genes, function(x){
    lapply(x,trimNAs)
})


list[geneDat, exp] = n_expressoExpr %>% sepExpr
design = n_expressoSamples

dpylrFriendly = cbind(design, t(exp))

dpylrFriendly %<>% arrange(MajorType,Neurotransmitter,ShinyNames) %>% filter(!is.na(PyramidalDeep))

genes %<>% lapply(function(x){
    names(x) = publishableNameDictionary$ShinyNames[match(names(x),publishableNameDictionary$PyramidalDeep)]
    return(x)
})


regionSamples = memoReg(design = dpylrFriendly,regionNames = 'Region',groupNames = 'ShinyNames',regionHierarchy=regionHierarchy)

names(regionSamples) = gsub('_.*','',names(regionSamples))
list[design,exp] = dpylrFriendly %>% sepExpr

colors = cellColors()

dir.create('analysis/01.SelectGenes/GenePlotsTop')
topGenes = genes[names(regionSamples)]



for (i in 1:len(regionSamples)){
    genes = topGenes[[i]]
    genes = genes[regionSamples[[i]] %>% unique %>% trimNAs]
    
    cellTypes = regionSamples[[i]] %>%trimNAs %>%  unique()
    relExp = exp[!is.na(regionSamples[[i]]),match(unlist(genes), geneDat$Gene.Symbol)]
    relGene = geneDat[match(unlist(genes), geneDat$Gene.Symbol),]
    relExp = apply(relExp,2,scale01)
    
    
    geneCellTypes = str_extract(names(unlist(genes)) , regexMerge(cellTypes,exact=TRUE))
    
    heatCols = toColor(design$ShinyNames[!is.na(regionSamples[[i]])], colors)
    geneCols = toColor(geneCellTypes, colors)
    
    png(paste0('analysis/01.SelectGenes/GenePlotsTop/',names(regionSamples)[i],'.png'),height=1400,width=1600)
    # heatmap.2(t(relExp),Rowv=F,Colv=F,trace='none',col=viridis(20),
    #           ColSideColors=heatCols$cols,RowSideColors=geneCols$cols,labRow='',labCol='', main = names(regionSamples[i]))
    # legend(title= 'Cell Types','bottomleft' ,legend= names(heatCols$palette), fill = heatCols$palette )
    
    anotCol = data.frame(Samples = design$ShinyNames[!is.na(regionSamples[[i]])])
    rownames(anotCol) = colnames(t(relExp))
    
    anotRow = data.frame('Specific Genes' = geneCellTypes, check.names=F)
    colnames(relExp) = genes %>% unlist 
    rownames(anotRow) = colnames((relExp))
    
    pheatmap(t(relExp),fontsize=30,gaps_row= anotRow$`Specific Genes` %>% duplicated %>% not %>% which %>% {.-1},
             show_rownames=TRUE ,
             show_colnames=FALSE,
             annotation_col=anotCol,
             annotation_row = anotRow,
             annotation_colors=list(Samples = heatCols$palette,
                                    'Specific Genes' = geneCols$palette),
             annotation_legend= T,
             cluster_rows=F,
             cluster_cols=F,
             main=paste(names(regionSamples[i]), 'marker gene expression'),
             color=viridis(20))
    dev.off()
}