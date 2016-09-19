library(ogbox)
library(oligoClasses)
library(affy)
devtools::load_all()
# download the parkinsons disease dasaets

# Lesnick et al -------------------------------------------
gseDown('GSE7621',outDir='data-raw/cel/GPL570/')
ogbox::getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL570',annotType='noParents')

dir.create('data-raw/LesnickParkinsons', showWarnings=FALSE)
softDown('GSE7621','data-raw/LesnickParkinsons/GSE7621_family.soft.gz')
system('gunzip data-raw/LesnickParkinsons/GSE7621_family.soft.gz')
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

cels = paste0('data-raw/cel/GPL570/', softData$GSM, '.cel')

affy = ReadAffy(filenames = cels)

norm = rma(affy)
annotated = gemmaAnnot(norm, 'data-raw/GemmaAnnots/GPL570')
annotated = mostVariable(annotated,threshold = 0)
write.csv(annotated, 'data-raw/LesnickParkinsons/GSE7621_parkinsonsExp.csv', row.names = F)


LesnickParkinsonsMeta = read.design('data-raw/LesnickParkinsons/GSE7621_parkinsonsMeta.tsv')
LesnickParkinsonsExp = read.exp('data-raw/LesnickParkinsons/GSE7621_parkinsonsExp.csv')

devtools::use_data(LesnickParkinsonsMeta, overwrite=TRUE)
devtools::use_data(LesnickParkinsonsExp, overwrite=TRUE)


# Moran et al dataset -------------------
gseDown('GSE8397',regex="A chip",outDir='data-raw/cel/GPL96/')
gseDown('GSE8397',regex="B chip",outDir='data-raw/cel/GPL97/')
ogbox::getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL96',annotType='noParents')
ogbox::getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL97',annotType='noParents')

dir.create('data-raw/MoranParkinsons', showWarnings=FALSE)
softDown('GSE8397',file='data-raw/MoranParkinsons/GSE8397_family.soft.gz')
system('gunzip data-raw/MoranParkinsons/GSE8397_family.soft.gz')
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
expTable = cbind(allGenes,allExp) %>% mostVariable(threshold=0)
write.csv(expTable, paste0('data-raw/MoranParkinsons/','GSE8397','_exp.csv'), row.names = F)
softData %>%
    filter(Platform=='GPL96') %>%
    select(-GSM) %>%
    write.design(file= paste0('data-raw/MoranParkinsons/','GSE8397','_des.tsv'))

MoranParkinsonsMeta = read.design('data-raw/MoranParkinsons/GSE8397_des.tsv')
MoranParkinsonsExp = read.exp('data-raw/MoranParkinsons/GSE8397_exp.csv')

devtools::use_data(MoranParkinsonsMeta, overwrite=TRUE)
devtools::use_data(MoranParkinsonsExp,overwrite=TRUE)
