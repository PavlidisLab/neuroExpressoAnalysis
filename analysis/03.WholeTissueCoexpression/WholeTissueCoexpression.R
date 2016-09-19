library(ogbox)
library(data.table)
library(stringr)
library(corpcor)
library(magrittr)
devtools::load_all()
library(homologene)



# handcrafting sets for coexpression
expObject = list(all = 'trabzuniRegionsExp',
                 cortex = 'trabzuniRegionsExp',
                 substantiaNigra= 'trabzuniRegionsExp',
                 cerebellum = 'trabzuniRegionsExp',
                 hippocampus = 'trabzuniRegionsExp',
                 thalamus = 'trabzuniRegionsExp',
                 chenCortex = 'chenExpr',
                 stanley1Cortex = 'stanleyStud1',
                 stanley3Cortex = 'stanleyStud3',
                 stanleyStud5 = 'stanleyStud5',
                 stanleyStud7 = 'stanleyStud7')

sets = list(all = rep(T,nrow(trabzuniRegionsMeta)),
            cortex = trabzuniRegionsMeta$brainRegion %in% c('frontal cortex', 
                                                            'occipital cortex',
                                                            'temporal cortex'),
            substantiaNigra = trabzuniRegionsMeta$brainRegion %in% 'substantia nigra',
            cerebellum = trabzuniRegionsMeta$brainRegion %in% 'cerebellar cortex',
            hippocampus = trabzuniRegionsMeta$brainRegion %in% 'hippocampus',
            thalamus = trabzuniRegionsMeta$brainRegion %in% 'thalamus',
            chenCortex = rep(T,ncol(chenExpr)),
            stanley1Cortex = stanleyMeta1$Profile %in% 'Cont',
            stanley3Cortex = stanleyMeta3$Profile %in% 'Cont',
            stanley5Cortex = stanleyMeta5$Profile %in% 'Cont',
            stanley7Cortex = stanleyMeta7$Profile %in% 'Cont')


# get all genes from all regions
allGenes = mouseMarkerGenes
allGenes = lapply(allGenes,lapply,function(x){
    mouse2human(x)$humanGene
})


# handcrafting gene sets
geneSets = list(all = allGenes[[1]][c('Astrocyte', 'Microglia','Oligo')],
                #all = allGenes[[1]],
                cortex = allGenes[['Cortex']] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                substantiaNigra = allGenes$Midbrain %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                cerebellum = allGenes$Cerebellum %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                hippocampus = allGenes$Hippocampus %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                thalamus = allGenes$Thalamus %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                chenCortex = allGenes[['Cortex']] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex1 = allGenes[['Cortex']]%>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex3 = allGenes[['Cortex']] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex5 = allGenes[['Cortex']] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex7 = allGenes[['Cortex']] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]}
)

# actual analysis

dupResolve = T
ps = lapply(1:len(sets), function(i){
    print(i)
    # take the relevant expression profile data
    expression = teval(expObject[[i]])
    # subset the human expression data to take in only samples from specific regions  
    setExpr =  expression[,sets[[i]]]
    rownames(setExpr) = rownames(expression)

    
    genes = geneSets[[i]]
    genes %<>% lapply(function(geneSub){geneSub = geneSub[geneSub %in% rn(setExpr)]})
    #expression[genes$GabaPV,] %>%t %>% cor
    #pearCor = setExpr[geneSets[[i]]$GabaPV,] %>% t %>% cor(method='spearman')
    # pearCor[geneSets[[i]]$GabaPV,geneSets[[i]]$GabaPV] %>% sm2vec %>% density %>% plot   
    #pearCor[tempGenes,tempGenes] %>% sm2vec %>% density %>% plot   
    medianExps = setExpr %>% apply(1,median)
    simuGenes = genes %>% lapply(function(geneSub){sapply(geneSub,selectRandom,500, medianExps)})
    
    # resolve duplicates
    if (dupResolve==T){
        simuGenes = lapply(1:len(simuGenes),function(j){
            simuGenes = simuGenes[[j]]
            if(len(simuGenes) == 0){
                return(list())
            } else if(ncol(simuGenes) == 1){
                return(simuGenes)
            }
            while (any(apply(apply(simuGenes,1,duplicated),2,any))){
                print(paste("had to resolve equality",names(genes)[j],'in',
                            sum(apply(apply(simuGenes,1,duplicated),2,any))))
                simuGenes[apply(apply(simuGenes,1,duplicated),2,any),] = 
                    t( apply(simuGenes[apply(apply(simuGenes,1,duplicated),2,any),,drop=F],1, function(x){
                        x = sample(x,len(x), replace = F)
                        x[duplicated(x)] = sapply(1:sum(duplicated(x)), function(y){
                            selectRandom(x[duplicated(x)][y],n=1, criteriaValue = medianExps,
                                         invalids=x[!x %in% x[duplicated(x)][y]])
                        })
                        return(x)
                    }))
            }
            return(simuGenes)
        })
        names(simuGenes) = names(genes)
    }
    
    # create the corr matrix from only the necesarry genes. hopefully will be slightly faster (it indeed is faster)
    realCors = genes %>% lapply(function(gene){
        tempCor = setExpr[gene,] %>% t %>% cor %>% sm2vec #%>% median
        #corMat[gene,gene] %>% sm2vec
    })
    
    simuCoexp = simuGenes %>% lapply(function(simuGene){
        if(len(simuGene)==0){
            return(NA)
        }
        simuGene %>% apply(1, function(x){
            tempCor = setExpr[x,] %>% t %>% cor %>% sm2vec # %>% median
        })
    })
    
    ps = sapply(1:len(simuCoexp), function(j){
        print(j)
        if(len(realCors[[j]])==0){
            return(NA)
        }
        print('dunnit')
        wilcox.test(simuCoexp[[j]] %>% as.vector, realCors[[j]], alternative = 'less')$p.value
    })
    
    ps[sapply(genes,len)<2] = NA
    
    #     
    #     ps = sapply(1:len(realCors), function(j){
    #         1-ecdf(simuCoexp[[j]])(realCors[j])
    #     })
    names(ps) = names(realCors)
    return(data.frame(ps,geneCounts = genes %>% sapply(len)))
})

names(ps) = names(geneSets)

pAdjusted  = relist(flesh=p.adjust(ps %>% lapply(function(x){x$ps}) %>% unlist,
                                   method='fdr'),
                    skeleton=ps %>% lapply(function(x){x$ps}))
ps = mapply(function(ps,pAdjusted){
    ps$ps = pAdjusted
    return(ps)},
    ps,pAdjusted,SIMPLIFY =FALSE)

dir.create('analysis//03.WholeTissueCoexpression/pValues', showWarnings=FALSE)
lapply(1:len(ps),function(i){    
    write.table(ps[[i]] %>% as.data.frame,col.names=F, sep= '\t',quote=F,file=paste0('analysis//03.WholeTissueCoexpression//pValues/',names(ps)[i]))
})

# table generation ------------
stanleyCortex1 = fread('analysis//03.WholeTissueCoexpression//pValues/stanleyCortex1',data.table=F)
stanleyCortex3 = fread('analysis//03.WholeTissueCoexpression//pValues/stanleyCortex3',data.table=F)
stanleyCortex5 = fread('analysis//03.WholeTissueCoexpression//pValues/stanleyCortex5',data.table=F)
stanleyCortex7 = fread('analysis//03.WholeTissueCoexpression//pValues/stanleyCortex7',data.table=F)

chenCortex = fread('analysis//03.WholeTissueCoexpression//pValues/chenCortex',data.table=F)
Trabzuni = fread('analysis//03.WholeTissueCoexpression//pValues/cortex',data.table=F)

cortex = data.frame('Cell Types'= stanleyCortex1$V1, 
                    'Trabzuni Dataset' = Trabzuni$V2,
                    'Trabzuni genecount' = Trabzuni$V3,
                    'Chen Dataset' = chenCortex$V2,
                    'Chen genecount' = chenCortex$V3,
                    'Stanley - AltarA' = stanleyCortex1$V2,
                    'AltarA genecount' =stanleyCortex1$V3,
                    'Stanley - Bahn' = stanleyCortex3$V2,
                    'Bahn genecount' =stanleyCortex3$V3,
                    'Stanley - Dobrin' = stanleyCortex5$V2,
                    'Dobrin genecount' =stanleyCortex5$V3,
                    'Stanley - Kato' = stanleyCortex7$V2,
                    'Kato genecount' =stanleyCortex7$V3,
                    check.names=F)

dir.create('analysis//03.WholeTissueCoexpression/publishTable', showWarnings=FALSE)
cortex[,-1] %<>% fixTable
cortex[,1] = publishableNameDictionary[match(cortex[,1],publishableNameDictionary$PyramidalDeep),]$ShinyNames
cortex %>% write.design('analysis//03.WholeTissueCoexpression/publishTable/cortexTable.tsv')

UCLCerebellum = fread('analysis//03.WholeTissueCoexpression/pValues/cerebellum',data.table=F)
UCLCerebellum[-1] %<>% fixTable
UCLCerebellum[,1] = publishableNameDictionary[match(UCLCerebellum[,1],publishableNameDictionary$PyramidalDeep),]$ShinyNames

UCLHippocampus = fread('analysis//03.WholeTissueCoexpression/pValues/hippocampus',data.table=F)
UCLHippocampus[-1] %<>% fixTable
UCLHippocampus[,1] = publishableNameDictionary[match(UCLHippocampus[,1],publishableNameDictionary$PyramidalDeep),]$ShinyNames

UCLSubstantia = fread('analysis//03.WholeTissueCoexpression/pValues/substantiaNigra',data.table=F)
UCLSubstantia[-1] %<>% fixTable
UCLSubstantia[,1] = publishableNameDictionary[match(UCLSubstantia[,1],publishableNameDictionary$PyramidalDeep),]$ShinyNames

UCLall = fread('analysis//03.WholeTissueCoexpression/pValues/all',data.table=F)
UCLall[-1] %<>% fixTable
UCLall[,1] = publishableNameDictionary[match(UCLall[,1],publishableNameDictionary$PyramidalDeep),]$ShinyNames

UCLThalamus = fread('analysis//03.WholeTissueCoexpression/pValues/thalamus',data.table=F)
UCLThalamus[-1] %<>% fixTable
UCLThalamus[,1] = publishableNameDictionary[match(UCLThalamus[,1],publishableNameDictionary$PyramidalDeep),]$ShinyNames

UCLCerebellum %>% write.design('analysis//03.WholeTissueCoexpression/publishTable/cerebellumTable.tsv')
UCLHippocampus %>% write.design('analysis//03.WholeTissueCoexpression/publishTable/HippocampusTable.tsv')
UCLSubstantia %>% write.design('analysis//03.WholeTissueCoexpression/publishTable/snTable.tsv')
UCLall %>% write.design('analysis//03.WholeTissueCoexpression/publishTable/allTable.tsv')
UCLThalamus %>% write.design('analysis//03.WholeTissueCoexpression/publishTable/thalamusTable.tsv')


