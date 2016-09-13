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

if(length(commandArgs(trailingOnly=TRUE))==0){
    start = 1
    end = 500
    secondChip = TRUE
} else{
    args <- commandArgs(trailingOnly = TRUE)
    start = as.numeric(args[1])
    end = as.numeric(args[2])
    secondChip = as.logical(args[3])
}


if (start == 1){
    # this is a quick way to select "goog enough" markers without doing permutations
    # output of this will not be robust to outliers. These genes are not used in the study
    # and are not readily available in the package
    markerCandidates(design = n_expressoSamples,
                     expression = n_expressoExpr,
                     outLoc = 'analysis//01.SelectGenes/Quick',
                     groupNames = 'PyramidalDeep',
                     regionNames = 'Region',
                     cores=16,
                     regionHierarchy = regionHierarchy
    )
    
    # quickly select genes exlusively for the samples from GPL1261. These genes are not used in the study and are not
    # readily available in the package
    if (secondChip){
        markerCandidates(design = 
                             ogbox::read.design('data-raw/Mouse_Cell_Type_Data/n_expressoSamples2.tsv'),
                         expression =
                             ogbox::read.exp('data-raw/Mouse_Cell_Type_Data/n_expressoExpr2.csv'),
                         outLoc = 'analysis//01.SelectGenes/Quick2',
                         groupNames = 'PyramidalDeep',
                         regionNames = 'Region',
                         cores=16,
                         regionHierarchy = regionHierarchy
        )
    }
}

# here we do the permutations required for selection of marker genes robust to outliers.
if(start==1){
    file.create('analysis//01.SelectGenes/Rotation/progress')
}
for (i in start:end){
    print(i)
    markerCandidates(design = n_expressoSamples,
                     expression = n_expressoExpr,
                     outLoc = paste0('analysis//01.SelectGenes/Rotation/',i),
                     groupNames = 'PyramidalDeep',
                     regionNames = 'Region',
                     rotate=0.33,
                     regionHierarchy = regionHierarchy,
                     cores=16
    )
    
}
cat(paste(start,end,'\n'),file='analysis//01.SelectGenes/Rotation/progress',append=TRUE)

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
                         groupNames = 'PyramidalDeep',
                         regionNames = 'Region',
                         rotate=0.33,
                         cores=16,
                         regionHierarchy = regionHierarchy
        )
    }
    cat(paste(start,end,'\n'),file='analysis//01.SelectGenes/Rotation2/progress',append=TRUE)
    
}

# if this is the last rotation, calculate the selection percentages of genes.
if(end == 500){
    # wait for all other branches to complete operation
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
                 cores = 16,foldChange = 0)
    if(secondChip){
        rotateSelect(rotationOut='analysis//01.SelectGenes/Rotation/',
                     rotSelOut='analysis/01.SelectGenes/RotSel2',
                     cores = 16, foldChange = 0)
    }
    
    
    # upon calculation of selection percentages in permutations, create a directory that houses genes
    # that are selected in more than 95% of the permutations
    allGenes = list(genes1 = pickMarkersAll('analysis/01.SelectGenes/RotSel/'))
    if (secondChip){
        allGenes = list(genes1 = allGenes[[1]],
                        genes2 = pickMarkersAll('analysis/01.SelectGenes/RotSel/'))
    }
    
    for (n in 1:len(allGenes)){
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
    # number of genes removed from microglia is needed in the paper
    genes = pickMarkersAll('analysis//01.SelectGenes/FinalGenes1/PyramidalDeep/')
    allMicroglia = genes %>% lapply(function(x){
        x['Microglia']
    }) %>% unlist %>% unique %>% len
    print(paste0('Microglia used to have ', allMicroglia, ' genes'))
    
    microglialException('analysis/01.SelectGenes/FinalGenes1/',cores=8)
    if (secondChip){
        microglialException('analysis/01.SelectGenes/FinalGenes2/',cores=8)
    }
    
    genes = pickMarkersAll('analysis//01.SelectGenes/FinalGenes1/PyramidalDeep/')
    allMicroglia = genes %>% lapply(function(x){
        x['Microglia']
    }) %>% unlist %>% unique %>% len
    print(paste0('Microglia now have ', allMicroglia, ' genes'))
    
    
    # after everything is done save the genes to the package
    mouseMarkerGenes = pickMarkersAll('analysis/01.SelectGenes/FinalGenes1/PyramidalDeep/')
    # Lpl is manually removed from the list as it is known to be expressed in adipocytes yet are not present
    # in our dataset
    mouseMarkerGenes %<>% lapply(function(x){
        x %>% lapply(function(y){
            y[!y %in% c('Lpl')]
        })
    })
    devtools::use_data(mouseMarkerGenes, overwrite=TRUE)
    
}
