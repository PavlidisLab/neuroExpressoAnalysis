library(ogbox)
library(oligoClasses)
library(affy)
devtools::load_all()
# download the parkinsons disease dasaets

downloadData = FALSE
downloadGemma = TRUE

# Lesnick et al -------------------------------------------
if (downloadGemma){
    ogbox::getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL570',annotType='noParents', overwrite =TRUE)
}

dir.create('data-raw/LesnickParkinsons', showWarnings=FALSE)
dir.create('data-raw/cel/GPL570/', showWarnings=FALSE)

if(downloadData){
    gseDown('GSE7621',outDir='data-raw/cel/GPL570/')
    softDown('GSE7621','data-raw/LesnickParkinsons/GSE7621_family.soft.gz')
    system('gunzip data-raw/LesnickParkinsons/GSE7621_family.soft.gz')
}
softData = softParser(softFile='data-raw/LesnickParkinsons/GSE7621_family.soft',expression=F)
softData = softData[c('!Sample_characteristics_ch1',
                      '!Sample_geo_accession',
                      '!Sample_title')]

names(softData) = c('Characteristic',
                    'GSM',
                    'title')

softData$scanDate = sapply(softData$GSM, function(x){
    celfileDate(paste0('data-raw/cel/GPL570/',x, '.cel'))
})

softData %<>% mutate(parkinson = grepl('PD',title)) %>% 
    mutate(female = grepl('female',Characteristic))


write.design(softData,'data-raw/LesnickParkinsons/GSE7621_parkinsonsMeta.tsv')
softData = read.design('data-raw/LesnickParkinsons/GSE7621_parkinsonsMeta.tsv')
cels = paste0('data-raw/cel/GPL570/', softData$GSM, '.cel')

affy = ReadAffy(filenames = cels)

norm = affy::rma(affy)
annotated = gemmaAnnot(norm, 'data-raw/GemmaAnnots/GPL570')
# save(annotated,file = 'data/lesnickPreLowExpression.rda')

# gender fix. nothing to fix here but leaving the code here just in case
bioGen = bioGender(annotated) %>% replaceElement(c('M' = FALSE,'F' = TRUE)) %$% newVector
annotated = annotated[,!cn(annotated) %in% (names(which(bioGen != softData$female))) ]
softData = softData[bioGen == softData$female,]

medExp = annotated %>% sepExpr %>% {.[[2]]} %>% unlist %>% median
annotated = mostVariable(annotated,threshold = medExp, threshFun= median)
write.csv(annotated, 'data-raw/LesnickParkinsons/GSE7621_parkinsonsExp.csv', row.names = F)


LesnickParkinsonsMeta = read.design('data-raw/LesnickParkinsons/GSE7621_parkinsonsMeta.tsv')
LesnickParkinsonsExp = read.exp('data-raw/LesnickParkinsons/GSE7621_parkinsonsExp.csv')

devtools::use_data(LesnickParkinsonsMeta, overwrite=TRUE)
devtools::use_data(LesnickParkinsonsExp, overwrite=TRUE)


# Moran et al dataset -------------------
dir.create('data-raw/cel/GPL97', showWarnings=FALSE)
dir.create('data-raw/cel/GPL96', showWarnings=FALSE)
dir.create('data-raw/MoranParkinsons', showWarnings=FALSE)

if(downloadData){
    gseDown('GSE8397',regex="A chip",outDir='data-raw/cel/GPL96/')
    gseDown('GSE8397',regex="B chip",outDir='data-raw/cel/GPL97/')
    softDown('GSE8397',file='data-raw/MoranParkinsons/GSE8397_family.soft.gz')
    system('gunzip data-raw/MoranParkinsons/GSE8397_family.soft.gz')
    
    
}
if(downloadGemma){
    ogbox::getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL96',annotType='noParents', overwrite = TRUE)
    ogbox::getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL97',annotType='noParents', overwrite = TRUE)
}
softData = softParser(softFile='data-raw/MoranParkinsons/GSE8397_family.soft',expression=F)
softData = softData[c("!Sample_characteristics_ch1 = age",
                      "!Sample_title",
                      "!Sample_source_name_ch1",
                      '!Sample_platform_id')]
softData$GSM = rownames(softData)
names(softData) = c('Age',
                    'Title',
                    'Source',
                    'Platform',
                    'GSM')

softData %<>% mutate(Sex = Age %>% 
                         as.char %>% 
                         strsplit(':') %>% 
                         sapply(function(x){x[2]}),
                     Age = Age %>%
                         as.char %>% 
                         strsplit(':') %>% 
                         sapply(function(x){x[1]}) %>%
                         gsub('; gender','',x=.) %>% as.numeric,
                     Title =Title %>% str_extract('.*?(?= -(| )[AB] chip)'),
                     Region = Title %>% str_extract(regexMerge(c('Superior frontal gyrus','Lateral substantia nigra','Medial substantia nigra'))),
                     Disease = Title %>% grepl('control',x=.) %>% toColor(palette=c('TRUE' = 'control','FALSE' = 'PD')) %$% cols,
                     patient = paste0(Disease, str_extract(Title,pattern='[0-9]+')))



readA = softData %>% filter(Platform=='GPL96') %>% select(GSM) %>% unlist
readB = softData %>% filter(Platform=='GPL97') %>% select(GSM) %>% unlist
affyA = ReadAffy(filenames=paste0('data-raw//cel//GPL96/',readA,'.cel' ))
affyB = ReadAffy(filenames=paste0('data-raw//cel//GPL97/',readB,'.cel' ))
affyA = affy::rma(affyA)
affyB = affy::rma(affyB)
affyA = gemmaAnnot(affyA, 'data-raw/GemmaAnnots/GPL96')
affyB = gemmaAnnot(affyB, 'data-raw/GemmaAnnots/GPL97')


aName = softData %>% filter(Platform=='GPL96') %>% select(Title) %>% unlist 
bName = softData %>% filter(Platform=='GPL97') %>% select(Title) %>% unlist
list[aGene,aExp] = affyA %>% sepExpr
list[bGene,bExp] = affyB %>% sepExpr
names(aExp) = aName
names(bExp) = bName
bExp = bExp[, match(aName,bName)]
allGenes = rbind(aGene,bGene)
allExp = rbind(aExp,bExp)

# nothing to fix here too
bioGen = bioGender(cbind(allGenes,allExp))
allExp = allExp[,!cn(allExp) %in% (names(which(bioGen != softData$Sex))) ]
softData = softData[bioGen == softData$Sex,]


medExp = allExp %>% unlist %>% median
expTable = cbind(allGenes,allExp) %>% mostVariable(threshold=medExp,
                                                   threshFun=median)
write.csv(expTable, paste0('data-raw/MoranParkinsons/','GSE8397','_exp.csv'), row.names = F)
softData %>%
    filter(Platform=='GPL96') %>%
    select(-GSM) %>%
    write.design(file= paste0('data-raw/MoranParkinsons/','GSE8397','_des.tsv'))

MoranParkinsonsMeta = read.design('data-raw/MoranParkinsons/GSE8397_des.tsv')
MoranParkinsonsExp = read.exp('data-raw/MoranParkinsons/GSE8397_exp.csv')

devtools::use_data(MoranParkinsonsMeta, overwrite=TRUE)
devtools::use_data(MoranParkinsonsExp,overwrite=TRUE)

# Zhang et al. dataset ----------------------------------------------
dir.create('data-raw/ZhangParkinsons', showWarnings=FALSE)
dir.create('data-raw/cel/GPL96', showWarnings=FALSE)

if(downloadData){
    gseDown('GSE20295',outDir='data-raw/cel/GPL96/')
    softDown('GSE20295',file='data-raw/ZhangParkinsons/GSE20295_family.soft.gz')
    system('gunzip data-raw/ZhangParkinsons/GSE20295_family.soft.gz')
}

if(downloadGemma){
    ogbox::getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL96',annotType='noParents', overwrite=TRUE)
}
softData = softParser(softFile='data-raw/ZhangParkinsons/GSE20295_family.soft',expression=F)
softData = softData[,c('!Sample_characteristics_ch1 = age',
                       '!Sample_characteristics_ch1 = brain region',
                       "!Sample_characteristics_ch1 = disease state",
                       "!Sample_characteristics_ch1 = gender",
                       '!Sample_geo_accession',
                       "!Sample_platform_id",
                       "!Sample_title")]
names(softData) = c('Age',
                    'brainRegion',
                    'diseaseState',
                    'gender',
                    'GSM',
                    'platform',
                    'title')

softData$patient = str_extract(softData$title, '^(T|[0-9]).*?(?= )')
softData$scanDate = lapply(softData$GSM, function(x){
    celfileDate(paste0('data-raw/cel/GPL96/',x, '.cel')) %>% as.Date(format = '%Y-%m-%d')
}) %>% unlist
cels = paste0('data-raw/cel/GPL96/', softData$GSM, '.cel')


affy = ReadAffy(filenames = cels)

norm = affy::rma(affy)
annotated = gemmaAnnot(norm, 'data-raw/GemmaAnnots/GPL96')
bioGender = bioGender(annotated)

softData$bioGender = bioGender %>% replaceElement(c("M" ='male','F' = 'female')) %$% newVector

annotated = annotated[cn(annotated) %>% str_extract('^.*?(?=\\.)') %>% 
                          {!. %in% softData$GSM[softData$bioGender != softData$gender]}]

softData = softData[softData$bioGender == softData$gender,]

write.design(softData,'data-raw/ZhangParkinsons/GSE20295_parkinsonsMeta.tsv')


medExp = annotated %>% sepExpr %>% {.[[2]]} %>% unlist %>% median

rowMedians = annotated %>% sepExpr %>% {.[[2]]} %>%
    apply(1, function(x){
        regions =  softData$brainRegion %>% unique
        sapply(regions, function(y){
            x[softData$brainRegion %in% y] %>% mean
        })
    }) %>% apply(2, max)

annotated =
    annotated[rowMedians > medExp,] %>%
    mostVariable(threshold=0)

write.csv(annotated, 'data-raw/ZhangParkinsons/GSE20295_parkinsonsExp.csv', row.names = F)


ZhangParkinsonsMeta = read.design('data-raw/ZhangParkinsons/GSE20295_parkinsonsMeta.tsv')
ZhangParkinsonsExp = read.exp('data-raw/ZhangParkinsons/GSE20295_parkinsonsExp.csv')

ZhangParkinsonsMeta$diseaseState %<>% replaceElement(c('Control' ='control',
                                                       "Parkinson's disease" = 'PD',
                                                       'Parkinsons disease' = "PD")) %$% newVector

list[gene, expression] = ZhangParkinsonsExp %>% sepExpr()

expressionControl = expression[ZhangParkinsonsMeta$diseaseState %in% 'control']
regionsControl = ZhangParkinsonsMeta$brainRegion[ZhangParkinsonsMeta$diseaseState %in% 'control']
outliersControl = sampleMixup(expressionControl,regionsControl)


expressionPD = expression[ZhangParkinsonsMeta$diseaseState %in% 'PD']
regionsPD = ZhangParkinsonsMeta$brainRegion[ZhangParkinsonsMeta$diseaseState %in% 'PD']
outliersPD = sampleMixup(expressionPD,regionsPD)

badSamples = c(gsub(pattern = '.cel',replacement = '',c(outliersControl,outliersPD)),
               ZhangParkinsonsMeta %>% filter(bioGender!=gender) %$% GSM)

ZhangParkinsonsExp = ZhangParkinsonsExp[!grepl(regexMerge(badSamples), colnames(ZhangParkinsonsExp))]
ZhangParkinsonsMeta %<>% filter(!GSM %in% badSamples)

devtools::use_data(ZhangParkinsonsMeta, overwrite=TRUE)
devtools::use_data(ZhangParkinsonsExp, overwrite=TRUE)
