devtools::load_all()

if(length(commandArgs(trailingOnly=TRUE))==0){
    start = 1
    end = 500
} else{
    args <- commandArgs(trailingOnly = TRUE)
    start = as.numeric(args[1])
    end = as.numeric(args[2])
}


if (start == 1){
    # this is a quick way to select "good enough" markers without doing permutations
    # output of this will not be robust to outliers. These genes are not used in the study
    # and are not readily available in the package
    markerCandidates(design = mouseBloodCellsSamples,
                     replicates='Repl',
                     sampleName='GSM',
                     expression =  mouseBloodCellsExp,
                     outLoc = 'analysis//06.BloodGeneProfiles/MouseMarkers/Quick',
                     groupNames = c('lm11','lm22'),
                     cores=2)
    
    markerCandidates(design = humanBloodCellsSamples,
                     sampleName ='sampleName',
                     expression =  humanBloodCellsExp,
                     outLoc = 'analysis//06.BloodGeneProfiles/HumanMarkers/Quick',
                     groupNames = c('lm11','lm22'),
                     cores=2)
  
}

# here we do the permutations required for selection of marker genes robust to outliers.
if(start==1){
    file.create('analysis//01.SelectGenes/Rotation/progress')
}
for (i in start:end){
    print(i)
    markerCandidates(design = mouseBloodCellsSamples,
                     replicates='Repl',
                     sampleName ='GSM',
                     expression =  mouseBloodCellsExp,
                     outLoc = paste0('analysis//06.BloodGeneProfiles/MouseMarkers/Rotation/',i),
                     groupNames = c('lm11','lm22'),
                     cores=2,
                     rotate = 0.33,
                     seed=i)
    
    markerCandidates(design = humanBloodCellsSamples,PMID='PubMed.ID',
                     sampleName ='sampleName',
                     expression =  humanBloodCellsExp,
                     outLoc = paste0('analysis//06.BloodGeneProfiles/HumanMarkers/Rotation/',i),
                     groupNames = c('lm11','lm22'),
                     cores=2,
                     rotate = 0.33,
                     seed=i)
}

cat(paste(start,end,'\n'),file='analysis//06.BloodGeneProfiles/progress',append=TRUE)


# if this is the last rotation, calculate the selection percentages of genes.
if(end == 500){
    # wait for all other branches to complete operation
    repeat{
        progress = read.table('analysis//06.BloodGeneProfiles/progress') %>% apply(1, function(x){
            x[1]:x[2]
        }) %>% sapply(len) %>% sum
        if (progress>=500){
            break
        }
        Sys.sleep(60) 
    }
    
    print('waiting complete')
    rotateSelect(rotationOut='analysis/06.BloodGeneProfiles/MouseMarkers/Rotation',
                 rotSelOut='analysis/06.BloodGeneProfiles/MouseMarkers/RotSel',
                 cores = 16,foldChange = 1)
    rotateSelect(rotationOut='analysis/06.BloodGeneProfiles/HumanMarkers/Rotation',
                 rotSelOut='analysis/06.BloodGeneProfiles/HumanMarkers/RotSel',
                 cores = 16, foldChange = 1)
    allGenes = list(MouseMarkers = pickMarkersAll('analysis/06.BloodGeneProfiles/MouseMarkers/RotSel'),
                    HumanMarkers = pickMarkersAll('analysis/06.BloodGeneProfiles/HumanMarkers/RotSel'))
    
    for (n in 1:len(allGenes)){
        genes = allGenes[[n]]
        for (i in 1:len(genes)){
            pieces = strsplit(names(genes)[i],'_')[[1]]
            if (is.na(pieces[2])){
                pieces[2] = pieces[1]
                pieces[1] ='All'
            }
            dir.create(paste0('analysis//06.BloodGeneProfiles/',names(allGenes)[n],'FinalGenes/',
                              pieces[2] , '/' , pieces[1]), 
                       showWarnings=F, recursive=T)
        }
    }
    mouseBloodMarkers = allGenes$MouseMarkers
    humanBloodMarkers = allGenes$HumanMarkers
    devtools::use_data(mouseBloodMarkers, overwrite=TRUE)
    devtools::use_data(humanBloodMarkers, overwrite=TRUE)
    
    commonGenesLm11= intersect(mouseBloodMarkers$lm11 %>% unlist %>% mouse2human %$% humanGene,
                               humanBloodMarkers$lm11 %>% unlist)
    commonGenesLm22= intersect(mouseBloodMarkers$lm22 %>% unlist %>% mouse2human %$% humanGene,
                               humanBloodMarkers$lm22 %>% unlist)
}
