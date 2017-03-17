# here we select marker genes provided in the object "mouseBrainMarkers" This script can run as a whole
# or can be run from the command line to divide up the rotation process.
# usage: Rscript geneSelection.R [start] [end] [secondChip]
# start and end are intigers that provide the interval of permutations
# 500 permutations are performed for this study
# if start = 1, quick selection will be performed by the script
# if end = 500 the script will wait for the output folder to be filled with 500
# permutations and finishes the gene selection process.
# second chip controls whether or not genes should be selected exclusively for GPL1261
# the default output location is data-raw. Other directories will be created if not already
# existing
devtools::load_all()
library(jsonlite)
library(stringi)
#library(markerGenesManuscript)

if(length(commandArgs(trailingOnly=TRUE))==0){
    start = 1
    end = 500
    firstChip = TRUE
    secondChip = FALSE
    signleCell = TRUE
} else{
    args <- commandArgs(trailingOnly = TRUE)
    start = as.numeric(args[1])
    end = as.numeric(args[2])
    firstChip = as.logical(args[5])
    if(is.na(firstChip)){
        firstChip = TRUE
    }
    secondChip = as.logical(args[3])
    singleCell = as.logical(args[4])
}

# processing single cell data ---------

# regionGroups = memoReg(n_expressoSamplesWithRNAseq,'Region',c('PyramidalDeep'),regionHierarchy)

# list[gene,exp] = n_expressoExprWithRNAseq %>% sepExpr
# cortexSamples = n_expressoSamplesWithRNAseq$sampleName[!is.na(regionGroups$Cortex_PyramidalDeep)]
# cortexMustSamples = n_expressoSamplesWithRNAseq$sampleName[!is.na(regionGroups$Cortex_PyramidalDeep) & 
#                                                                (!(n_expressoSamplesWithRNAseq$PyramidalDeep %in%  n_expressoSamples$PyramidalDeep & 
#                                                                       n_expressoSamplesWithRNAseq$Platform %in% 'RNAseq'))]
# 
# cortexOldSamples =  n_expressoSamplesWithRNAseq$sampleName[!grepl('Layer', n_expressoSamplesWithRNAseq$PyramidalDeep) & !is.na(regionGroups$Cortex_PyramidalDeep)]
# 
# n_expressoSamplesWithRNAseq$PyramidalDeepNoNewPyramidal = n_expressoSamplesWithRNAseq$PyramidalDeep
# n_expressoSamplesWithRNAseq$PyramidalDeepNoNewPyramidal[!n_expressoSamplesWithRNAseq$sampleName %in% cortexOldSamples] = NA
# 
# n_expressoSamplesWithRNAseq$PyramidalDeepNoSingleCellUnlessYouHaveTo =  n_expressoSamplesWithRNAseq$PyramidalDeep
# n_expressoSamplesWithRNAseq$PyramidalDeepNoSingleCellUnlessYouHaveTo[!n_expressoSamplesWithRNAseq$sampleName %in% cortexMustSamples] = NA
# 
# n_expressoSamplesWithRNAseq$PyramidalDeepNoSingleCellUnlessYouHaveToAndNoNewPyramidals = n_expressoSamplesWithRNAseq$PyramidalDeep
# n_expressoSamplesWithRNAseq$PyramidalDeepNoSingleCellUnlessYouHaveToAndNoNewPyramidals[!n_expressoSamplesWithRNAseq$sampleName %in% cortexMustSamples | !n_expressoSamplesWithRNAseq$sampleName %in% cortexOldSamples] = NA
# 
# meltedSingleCells$PyramidalDeepNoNewPyramidal = meltedSingleCells$PyramidalDeep
# meltedSingleCells$PyramidalDeepNoNewPyramidal[!meltedSingleCells$sampleName %in% cortexOldSamples] = NA


# quick selection ---------------------------
if (start == 1){
    # this is a quick way to select "goog enough" markers without doing permutations
    # output of this will not be robust to outliers. These genes are not used in the study
    # and are not readily available in the package
    if(singleCell){
        # markerCandidates(design = n_expressoSamplesWithRNAseq[n_expressoSamplesWithRNAseq$sampleName %in% cortexSamples,],
        #                  expression = cbind(gene,exp[cortexSamples]),
        #                  outLoc = 'analysis//01.SelectGenes/QuickSingleCell',
        #                  groupNames = c('PyramidalDeep',
        #                                 'PyramidalDeepNoNewPyramidal',
        #                                 'PyramidalDeepNoSingleCellUnlessYouHaveTo',
        #                                 'PyramidalDeepNoSingleCellUnlessYouHaveToAndNoNewPyramidals'),
        #                  #groupNames = 'DopaSelect',
        #                  #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
        #                  regionNames = NULL,
        #                  cores=16,
        #                  regionHierarchy= NULL)
        ptm <- proc.time()

        markerCandidates(design = meltedSingleCells,foldChangeThresh = 10,minimumExpression = 2.5,
                         expression = data.frame(Gene.Symbol = rn(TasicPrimaryMeanLog),TasicPrimaryMeanLog,check.names = FALSE),
                         outLoc = 'analysis//01.SelectGenes/QuickJustSingleCell',
                         groupNames = c('PyramidalDeep','CellTypes'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = NULL,
                         cores=16,
                         regionHierarchy= NULL)
        singleCellTime = proc.time() - ptm
        
        ptm <- proc.time()
        
        markerCandidates(design = meltedSingleCells,foldChangeThresh = 8,minimumExpression = 1.5,
                         expression = data.frame(Gene.Symbol = rn(TasicPrimaryMeanLog),TasicPrimaryMeanLog,check.names = FALSE),
                         outLoc = 'analysis//01.SelectGenes/QuickJustSingleCellRelax',
                         groupNames = c('PyramidalDeep','CellTypes'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = NULL,
                         cores=16,
                         regionHierarchy= NULL)
        singleCellTimeRelax = proc.time() - ptm
        
    }
    
    if(firstChip){
        ptm <- proc.time()
        markerCandidates(design = n_expressoSamples,
                         expression = n_expressoExpr,
                         outLoc = 'analysis//01.SelectGenes/Quick',
                         groupNames = c('PyramidalDeep','CellTypes'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = 'Region',
                         cores=16,
                         regionHierarchy = regionHierarchy)
        neuroExpTime = proc.time() - ptm
        
    }
    # quickly select genes exlusively for the samples from GPL1261. These genes are not used in the study and are not
    # readily available in the package
    if (secondChip){
        markerCandidates(design = 
                             ogbox::read.design('data-raw/Mouse_Cell_Type_Data/n_expressoSamples2.tsv'),
                         expression =
                             ogbox::read.exp('data-raw/Mouse_Cell_Type_Data/n_expressoExpr2.csv'),
                         outLoc = 'analysis//01.SelectGenes/Quick2',
                         groupNames = c('PyramidalDeep','BroadTypes','DopaSelect'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = 'Region',
                         cores=16,
                         regionHierarchy = regionHierarchy)
    }
}


# rotations ----------------------------
if(start==1){
    file.create('analysis//01.SelectGenes/Rotation/progress')
}
if(firstChip){
    for (i in start:end){
        print(i)
        markerCandidates(design = n_expressoSamples,
                         expression = n_expressoExpr,
                         outLoc = paste0('analysis//01.SelectGenes/Rotation/',i),
                         groupNames = c('PyramidalDeep','CellTypes'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = 'Region',
                         rotate=0.33,
                         regionHierarchy = regionHierarchy,
                         cores=16,
                         seed = i)
        
    }
}
cat(paste(start,end,'\n'),file='analysis//01.SelectGenes/Rotation/progress',append=TRUE)

# rotation with single cells ------------------------------
if(singleCell){
    if(start==1){
        file.create('analysis//01.SelectGenes/RotationSingleCell/progress')
    }
    for(i in start:end){
        # markerCandidates(design = n_expressoSamplesWithRNAseq[n_expressoSamplesWithRNAseq$sampleName %in% cortexSamples,],
        #                  expression = cbind(gene,exp[cortexSamples]),
        #                  outLoc = paste0('analysis//01.SelectGenes/RotationSingleCell/',i),
        #                  groupNames = c('PyramidalDeep',
        #                                 'PyramidalDeepNoNewPyramidal',
        #                                 'PyramidalDeepNoSingleCellUnlessYouHaveTo',
        #                                 'PyramidalDeepNoSingleCellUnlessYouHaveToAndNoNewPyramidals'),
        #                  #groupNames = 'DopaSelect',
        #                  #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
        #                  regionNames = NULL,
        #                  cores=16,
        #                  rotate=0.33,
        #                  regionHierarchy= NULL)
        
        markerCandidates(design = meltedSingleCells,foldChangeThresh = 10, minimumExpression = 2.5,
                         expression = data.frame(Gene.Symbol = rn(TasicPrimaryMeanLog),TasicPrimaryMeanLog,check.names = FALSE),
                         outLoc = paste0('analysis//01.SelectGenes/RotationJustSingleCell/',i),
                         groupNames = c('PyramidalDeep','CellTypes'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = NULL,
                         cores=15,
                         rotate = 0.33,
                         regionHierarchy= NULL)
        
        markerCandidates(design = meltedSingleCells,foldChangeThresh = 8,minimumExpression = 1.5,
                         expression = data.frame(Gene.Symbol = rn(TasicPrimaryMeanLog),TasicPrimaryMeanLog,check.names = FALSE),
                         outLoc = paste0('analysis//01.SelectGenes/RotationJustSingleCellRelaxed/',i),
                         groupNames = c('PyramidalDeep','CellTypes'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = NULL,
                         cores=4,
                         rotate = 0.33,
                         regionHierarchy= NULL)
        
        
    }
    cat(paste(start,end,'\n'),file='analysis//01.SelectGenes/RotationSingleCell/progress',append=TRUE)
}
# second chip rotations -------------------
if(secondChip){
    if(start==1){
        file.create('analysis//01.SelectGenes/Rotation2/progress')
    }
    for (i in start:end){
        print(i)
        markerCandidates(design = 
                             ogbox::read.design('data-raw/Mouse_Cell_Type_Data/n_expressoSamples2.tsv'),
                         expression = 
                             ogbox::read.exp('data-raw/Mouse_Cell_Type_Data/n_expressoExpr2.csv'),
                         outLoc = paste0('analysis//01.SelectGenes/Rotation2/',i),
                         groupNames = c('PyramidalDeep','BroadTypes','DopaSelect'),
                         #groupNames = 'DopaSelect',
                         #groupNames = c('AstroInactiveAlone','AstroReactiveAlone'),
                         regionNames = 'Region',
                         rotate=0.33,
                         cores=16,
                         regionHierarchy = regionHierarchy,
                         seed = i)
    }
    cat(paste(start,end,'\n'),file='analysis//01.SelectGenes/Rotation2/progress',append=TRUE)
    
}

# RotSel: if this is the last rotation, calculate the selection percentages of genes. ----------------
if(end == 500){
    # wait for all other branches to complete operation
    if(firstChip){
        repeat{
            progress = read.table('analysis//01.SelectGenes/Rotation/progress') %>% apply(1, function(x){
                x[1]:x[2]
            }) %>% sapply(len) %>% sum
            if (progress>=500){
                break
            }
            Sys.sleep(60) 
        }
        
        print('waiting complete')
        rotateSelect(rotationOut='analysis//01.SelectGenes/Rotation/',
                     rotSelOut='analysis/01.SelectGenes/RotSel',
                     cores = 16,foldChange = 1)
    }
    # rotsel second chip
    if(secondChip){
        repeat{
            progress = read.table('analysis//01.SelectGenes/Rotation2/progress') %>% apply(1, function(x){
                x[1]:x[2]
            }) %>% sapply(len) %>% sum
            if (progress>=500){
                break
            }
            Sys.sleep(60) 
        }
        rotateSelect(rotationOut='analysis//01.SelectGenes/Rotation/',
                     rotSelOut='analysis/01.SelectGenes/RotSel2',
                     cores = 16, foldChange = 1)
    }
    # rotsel single cells
    if(singleCell){
        repeat{
            progress = read.table('analysis//01.SelectGenes/RotationSingleCell//progress') %>% apply(1, function(x){
                x[1]:x[2]
            }) %>% sapply(len) %>% sum
            if (progress>=500){
                break
            }
            Sys.sleep(60) 
        }
        rotateSelect(rotationOut='analysis//01.SelectGenes/RotationSingleCell//',
                     rotSelOut='analysis/01.SelectGenes/RotSelSingleCell',
                     cores = 16, foldChange = 1)
        
        rotateSelect(rotationOut='analysis//01.SelectGenes/RotationJustSingleCell//',
                     rotSelOut='analysis/01.SelectGenes/RotSelJustSingleCell',
                     cores = 16, foldChange = 1)
        
        rotateSelect(rotationOut='analysis//01.SelectGenes/RotationJustSingleCellRelaxed//',
                     rotSelOut='analysis/01.SelectGenes/RotSelJustSingleCellRelaxed',
                     cores = 16, foldChange = 1)
    }
    # upon calculation of selection percentages in permutations, create a directory that houses genes -----
    # that are selected in more than 95% of the permutations
    if(firstChip){
        allGenes = list(genes1 = pickMarkersAll('analysis/01.SelectGenes/RotSel/'))
    } else
        allGenes$gene1 = NA
}
if (secondChip){
    allGenes = c(allGenes = allGenes[[1]],
                 list(genes3 = pickMarkersAll('analysis/01.SelectGenes/RotSel/')))
} else{
    allGenes$gene2=NA
}
if(singleCell){
    allGenes = c(allGenes, list(genes3 = pickMarkersAll('analysis/01.SelectGenes/RotSelSingleCell/')),
                 list(genes4 = pickMarkersAll('analysis/01.SelectGenes/RotSelJustSingleCell//')),
                 list(genes5 = pickMarkersAll(('analysis/01.SelectGenes/RotSelJustSingleCellRelaxed/'))))
}


for (n in 1:len(allGenes)){
    if(is.na(allGenes[[n]])){
        next
    }
    genes = allGenes[[n]]
    for (i in 1:len(genes)){
        pieces = strsplit(names(genes)[i],'_')[[1]]
        if (is.na(pieces[2])){
            pieces[2] = pieces[1]
            pieces[1] ='All'
        }
        dir.create(paste0('analysis//01.SelectGenes/FinalGenes',n,'/',
                          pieces[2] , '/' , pieces[1]), 
                   showWarnings=F, recursive=T)
        
        
        for (j in 1:len(genes[[i]])){
            write.table(genes[[i]][[j]],
                        paste0('analysis/01.SelectGenes/FinalGenes',n,'/',
                               pieces[2],'/',pieces[1],'/', 
                               names(genes[[i]])[j]),
                        row.names=F,
                        quote=F,
                        col.names=F      
            )
        }
    }
}
# here we do some wrangling of the gene list to deatl with astrocytes and microglia

# remove activated microglia genes

if(firstChip){
    # number of genes removed from microglia is needed in the paper
    genes = pickMarkersAll('analysis//01.SelectGenes/FinalGenes1/PyramidalDeep/')
    allMicroglia = genes %>% lapply(function(x){
        x['Microglia']
    }) %>% unlist %>% unique %>% len
    print(paste0('Microglia used to have ', allMicroglia, ' genes'))
    
    microglialException('analysis/01.SelectGenes/FinalGenes1/',cores=8)
    
    genes = pickMarkersAll('analysis//01.SelectGenes/FinalGenes1/PyramidalDeep/')
    allMicroglia = genes %>% lapply(function(x){
        x['Microglia']
    }) %>% unlist %>% unique %>% len
    print(paste0('Microglia now have ', allMicroglia, ' genes'))
}
if (secondChip){
    microglialException('analysis/01.SelectGenes/FinalGenes2/',cores=8)
}
if(singleCell){
    microglialException('analysis/01.SelectGenes/FinalGenes3/',cores=8)
    microglialException('analysis/01.SelectGenes/FinalGenes4/',cores=8)
    microglialException('analysis/01.SelectGenes/FinalGenes5/',cores=8)
}



if(firstChip){
    genes = pickMarkersAll('analysis//01.SelectGenes/FinalGenes1/PyramidalDeep/')
    
    allS100 = genes %>% lapply(function(x){
        x['Pyramidal_S100a10']
    }) %>% unlist %>% unique %>% len
    print(paste0('S100a10 pyramdials used to have ', allS100, ' genes'))
    
    s100a10exception('analysis/01.SelectGenes/FinalGenes1/',cores=8)
    genes = pickMarkersAll('analysis//01.SelectGenes/FinalGenes1/PyramidalDeep/')
    allS100 = genes %>% lapply(function(x){
        x['Pyramidal_S100a10']
    }) %>% unlist %>% unique %>% len
    print(paste0('S100a10 pyramdials now have ', allS100, ' genes'))
}
if (secondChip){
    s100a10exception('analysis/01.SelectGenes/FinalGenes2/',cores=8)
}
if(singleCell){
    # s100a10exception('analysis/01.SelectGenes/FinalGenes3/',cores=8)
    # s100a10exception('analysis/01.SelectGenes/FinalGenes4/',cores=8)
}


# Lpl is manually removed from the list as it is known to be expressed in adipocytes yet are not present
# in our dataset
bannedGenes = c('Lpl','S100a10')
if(firstChip){
    banGenes(restDir = 'analysis/01.SelectGenes/FinalGenes1/',
             bannedGenes= bannedGenes,
             cores=8)
}
if (secondChip){
    banGenes('analysis/01.SelectGenes/FinalGenes2/',bannedGenes = bannedGenes,cores=8)
}
if(singleCell){
    banGenes('analysis/01.SelectGenes/FinalGenes3/',bannedGenes = bannedGenes,cores=8)
    banGenes('analysis/01.SelectGenes/FinalGenes4/',bannedGenes = bannedGenes,cores=8)
    banGenes('analysis/01.SelectGenes/FinalGenes5/',bannedGenes = bannedGenes,cores=8)
}

if(firstChip){
    
    genes = pickMarkersAll('analysis//01.SelectGenes/FinalGenes1/PyramidalDeep/')
    
    assertthat::validate_that(all(!bannedGenes %in% unlist(genes)))
    
    mouseMarkerGenes = genes
    
    
    devtools::use_data(mouseMarkerGenes, overwrite=TRUE)
    
    toPlotGenes = mouseMarkerGenes %>% lapply(function(x){
        x = x[!grepl('Microglia_',names(x))]
        x %<>% lapply(function(y){
            y[!grepl('[|]', y)]
        })
    })
    
    
    toPlotGenes %<>% lapply(function(x){
        x %<>%sapply(len)
        x[cellOrder] %>% trimNAs
        x = x[!grepl('Microglia_',names(x))]
        names(x) = publishableNameDictionary$ShinyNames[match(names(x),publishableNameDictionary$PyramidalDeep)]
        return(x)
    })
    
    toPlotGenes$All = toPlotGenes$All[c('Astrocyte','Oligodendrocyte','Microglia')]
    toPlotGenes[-1] %<>% lapply(function(x){
        x = x[!names(x) %in% c('Astrocyte','Oligodendrocyte','Microglia')]
    })
    # take the bottom ones in the region tree
    rockBottom = regionHierarchy %>% unlist %>% names %>% str_extract(pattern='(?<=[.])([A-Z]|[a-z])*$')
    rockBottom = c(rockBottom,'Midbrain')
    toPlotGenes = toPlotGenes[c('All', rockBottom)]
    
    toPlotGenes %<>% lapply(function(x){
        sapply(1:len(x),function(i){
            genes = x[i]
            samples = n_expressoSamples %>% filter(ShinyNames %in% names(x[i])) %>% nrow
            sources = 
                n_expressoSamples %>% 
                filter(ShinyNames %in% names(x[i])) %>% 
                select(Reference,GSE) %>% 
                unique %>% {
                    paste0(.[,1],' (', .[,2], ')')
                } %>% paste(collapse = ', ')
            out = c(names(x[i]),samples, genes, sources )
            names(out) = NULL
            return(out)
        }) %>% as.data.frame %>% t
    }) 
    
    file.create('analysis/01.SelectGenes/geneTable.tsv')
    lapply(1:len(toPlotGenes), function(i){
        print(i)
        if(i == 1){
            append = FALSE
        } else {
            append = TRUE
        }
        cat(paste0(names(toPlotGenes)[i],'\n'),
            append = append,
            file= 'analysis/01.SelectGenes/geneTable.tsv')
        write.table(toPlotGenes[[i]], file = 'analysis/01.SelectGenes/geneTable.tsv',
                    sep = "\t",
                    quote = F, col.names = F,
                    row.names = F, append = TRUE)
    })
    
    # gene list in single files -------
    mouseMarkerGenes %>% toJSON(pretty=TRUE) %>% writeLines('analysis/01.SelectGenes/markerGenes.json')
    
    sheet = loadWorkbook('analysis/01.SelectGenes/markerGenes.xls', create = TRUE)
    
    dir.create('analysis/01.SelectGenes/markerGeneTSVs')
    1:len(mouseMarkerGenes) %>% sapply(function(i){
        out = stri_list2matrix(mouseMarkerGenes[[i]]) %>% as.data.frame
        names(out) = names(mouseMarkerGenes[[i]])
        write.table(out,file = file.path('analysis/01.SelectGenes/markerGeneTSVs',names(mouseMarkerGenes[i])),na= '', sep = "\t", quote = F, row.names = F)
        createSheet(sheet, name = names(mouseMarkerGenes[i]))
        writeWorksheet(sheet, out, sheet =  names(mouseMarkerGenes[i]), startRow = 1, startCol = 1)
    })
    saveWorkbook(sheet)
    
    
    # create archive
    system('rar -ep1 a analysis/01.SelectGenes/markerGenes.rar analysis/01.SelectGenes/FinalGenes1/PyramidalDeep/*')
}
