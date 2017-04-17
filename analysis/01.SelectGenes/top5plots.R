# assertthat::validate_that(all(!bannedGenes %in% unlist(genes)))
library(neuroExpressoAnalysis)
library(dplyr)
library(viridis)
library(reshape2)
library(ogbox)
library(scales)
library(gplots)
library(stringr)
library(magrittr)
library(pheatmap) # dev version is used
# expression of top 5 heatmap --------
order = cellOrder %>% translatePublishable()

genes = mouseMarkerGenes
#genes$Cortex %<>% c(list('Pyramidal' = mouseMarkerGenes$Cortex$Pyramidal))

genes %<>% lapply(function(x){x[!grepl(pattern='Microglia_',names(x))]})
geneNames = names(genes)

genesMicroarray = lapply(1:len(genes), function(i){
    nm = names(genes[[i]])
    nm = nm[nm %in% (n_expressoSamples %>% select(CellTypes,PyramidalDeep) %>% unlist %>% unique)]
    lapply(1:len(nm), function(j){
        print(paste(i,j))
        component1 = names(genes[i])
        component2 = names(genes[[i]][nm[j]])
        component1 = paste0(component1,'_')
        if (component1=='All_'){
            component1=''
        }
        if(component2 == 'Pyramidal'){
            scores = read.table(paste0('analysis/01.SelectGenes/Quick/', component1,'CellTypes/',component2))
        } else{
            scores = read.table(paste0('analysis/01.SelectGenes/Quick/', component1,'PyramidalDeep/',component2))
        }
        (scores %>% filter(V1 %in% genes[[i]][[nm[j]]]))[1:5,]$V1 %>% unlist %>% as.char
    }) ->out
    names(out) = nm
    return(out)
})

names(genesMicroarray) = geneNames
genesMicroarray = lapply(genesMicroarray, function(x){
    lapply(x,trimNAs)
})

genesSingleCell  = lapply(1:len(genes$Cortex), function(i){
    print(names(genes$Cortex[i]))
    if(names(genes$Cortex[i]) == 'Pyramidal'){
        scores = read.table(paste0('analysis/01.SelectGenes/QuickJustSingleCell//CellTypes/',names(genes$Cortex[i])))
    } else{
        scores = read.table(paste0('analysis/01.SelectGenes/QuickJustSingleCell//PyramidalDeep/',names(genes$Cortex[i])))
    }
    assertthat::are_equal(all(genes$Cortex[[i]] %in% scores$V1),TRUE)
    (scores %>% filter(V1 %in% genes$Cortex[[i]]))[1:5,]$V1 %>% unlist %>% as.char
})
names(genesSingleCell) = names(genes$Cortex)
genesSingleCell %<>% lapply(trimNAs)

n_expressoExpr %<>% filter(!grepl(pattern='\\|',Gene.Symbol))
list[geneDat, exp] = n_expressoExpr %>% sepExpr
rn(exp) =geneDat$Gene.Symbol
design = n_expressoSamples

dpylrFriendly = cbind(design, t(exp))

dpylrFriendly %<>%
    mutate(ShinyNames = CellTypes %>% translatePublishable %>% factor(levels = order[order %in% (CellTypes %>% translatePublishable)])) %>% 
    arrange(ShinyNames) %>% filter(!is.na(PyramidalDeep))

genesMicroarray %<>% lapply(function(x){
    names(x) = publishableNameDictionary$ShinyNames[match(names(x),publishableNameDictionary$PyramidalDeep)]
    return(x)
})
names(genesSingleCell) = publishableNameDictionary$ShinyNames[match(names(genesSingleCell),publishableNameDictionary$PyramidalDeep)]

regionSamples = memoReg(design = dpylrFriendly,regionNames = 'Region',groupNames = 'ShinyNames',regionHierarchy=regionHierarchy)

names(regionSamples) = gsub('_.*','',names(regionSamples))
list[design,exp] = dpylrFriendly %>% sepExpr
exp %<>% cbind(NA)

colors = cellColors()

dir.create('analysis/01.SelectGenes/GenePlotsTop', showWarnings = FALSE)
topGenes = genesMicroarray[names(regionSamples)]

null = list(NULL)
names(null) = NA

for (i in 1:len(regionSamples)){
    print(i)
    
    if(names(regionSamples)[i] == 'Cortex'){
        genes = names(genesMicroarray$Cortex) %>% sapply(function(x){
            c( genesMicroarray$Cortex[[x]],genesSingleCell[[x]]) %>% unique
        })
        
        genes = c(genes[order[order %in% c(design$ShinyNames %>% as.char,'Pyramidal')]] %>% trimElement(null),
                  genesSingleCell[order[order %in% names(genesSingleCell)[!names(genesSingleCell) %in% c(design$ShinyNames %>% as.char,'Pyramidal')]]])
        
    }  else{
        genes = topGenes[[i]] %>% {.[order[order %in% design$ShinyNames]]} %>% trimElement(null)
    }
    
    
   # genes = genes[order[order %in% names(genes)]]
    
    cellTypes = regionSamples[[i]] %>% as.char %>%trimNAs %>%  unique()
    if(names(regionSamples)[i] == 'Cortex'){
        cellTypes = c(cellTypes, meltedSingleCells$CellTypes %>% translatePublishable %>% as.char) %>% unique
    }
    relExp = exp[!is.na(regionSamples[[i]]),
                 match(unlist(genes), geneDat$Gene.Symbol) %>% 
                     replaceElement(c(a='b'),NAreplace = ncol(exp)) %$% 
                     newVector %>% 
                     as.numeric]
    
    relGene = geneDat[match(unlist(genes), geneDat$Gene.Symbol),]
    relExp = apply(relExp,2,scale01)
    cn(relExp) = unlist(genes)
    
    if(names(regionSamples)[i] == 'Cortex'){
        relExp2 = TasicPrimaryMeanLog[match(unlist(genes),rn(TasicPrimaryMeanLog)),] %>% t
        
        meltedSingleCells %<>%
            mutate(ShinyNames = CellTypes %>% translatePublishable %>% 
                       factor(
                           levels =  
                               c(order[order%in% (regionSamples$Cortex %>% unique %>% trimNAs)],
                                 order[order %in% (CellTypes %>% translatePublishable)[!(CellTypes %>% translatePublishable) %in% (regionSamples$Cortex %>% unique %>% trimNAs)]]))) %>% 
            arrange(ShinyNames)
        relExp2 = relExp2[meltedSingleCells$sampleName,] %>% apply(2,scale01)
        cn(relExp2) = unlist(genes)
        
        relExp = rbind(relExp2,relExp)
    }
    
    geneCellTypes = repIndiv(names(genes), genes %>% sapply(len))
    
    # geneCellTypes = str_extract(names(unlist(genes)) , regexMerge(order,exact=TRUE))
    
    heatCols = toColor(as.char(design$ShinyNames)[!is.na(regionSamples[[i]])], colors)
    if(names(regionSamples)[i] == 'Cortex'){
        heatCols = meltedSingleCells$ShinyNames %>% as.char %>% toColor(colors)
    }
    
    geneCols = toColor(geneCellTypes, colors)
    
    if(names(regionSamples)[i] == 'Cortex'){
        rexTasic = TasicPrimaryMeanLog[,]
    }
    
    
   png(paste0('analysis/01.SelectGenes/GenePlotsTop/',names(regionSamples)[i],'.png'),height=1900,width=2000)
    # heatmap.2(t(relExp),Rowv=F,Colv=F,trace='none',col=viridis(20),
    #           ColSideColors=heatCols$cols,RowSideColors=geneCols$cols,labRow='',labCol='', main = names(regionSamples[i]))
    # legend(title= 'Cell Types','bottomleft' ,legend= names(heatCols$palette), fill = heatCols$palette )
    
    anotCol = data.frame(Samples = design$ShinyNames[!is.na(regionSamples[[i]])])
    if(names(regionSamples)[i] == 'Cortex'){
        anotCol = rbind(data.frame(Samples= meltedSingleCells$ShinyNames),anotCol)
    }
    
    rownames(anotCol) = colnames(t(relExp))
    
    anotRow = data.frame('Specific Genes' = geneCellTypes, check.names=F)
    colnames(relExp) = genes %>% unlist 
    rownames(anotRow) = colnames((relExp))
    
    pheatmap(t(relExp),
             fontsize=30,
             #fontsize = 11,
             #gaps_col = sum(!is.na(regionSamples[[i]])),
             gaps_col = relExp2 %>% ncol,
             gaps_row= anotRow$`Specific Genes` %>% duplicated %>% not %>% which %>% {.-1},
             show_rownames=TRUE ,
             show_colnames=FALSE,
             annotation_col=anotCol,
             annotation_row = anotRow,
             annotation_colors=list(Samples = heatCols$palette[match(anotCol$Samples %>% unique,names(heatCols$palette))],
                                    'Specific Genes' = geneCols$palette[match(anotRow$`Specific Genes` %>% unique,names(heatCols$palette))]),
             annotation_legend= T,
             cluster_rows=F,
             cluster_cols=F,
             main=paste(names(regionSamples[i]), 'marker gene expression'),
             color=viridis(20),
             na_col = 'grey',
             border_color = NA)
    dev.off()
}

# single cell plot






