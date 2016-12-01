library(limma)
library(data.table)
library(stringr)
library(sva)
library(oligoClasses)
library(magrittr)
library(ogbox)
devtools::load_all()

# trabzuni dataset ----------------------
#ogbox::getGemmaAnnot('GPL5175','data-raw/GemmaAnnots/GPL5175',annotType='noParents')

dir.create('data-raw/TrabzuniRegions/',showWarnings=FALSE)
# download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60862/soft/GSE60862_family.soft.gz',
#               destfile='data-raw/TrabzuniRegions/GSE60862_family.soft.gz')
# system('gunzip data-raw/TrabzuniRegions/GSE60862_family.soft.gz')

# deal with softfile
softData = softParser('data-raw/TrabzuniRegions/GSE60862_family.soft')

###
softData = softData[,c('!Sample_characteristics_ch1 = age at death (in years)', 
                       '!Sample_characteristics_ch1 = ancestry',
                       '!Sample_characteristics_ch1 = brain bank',
                       '!Sample_characteristics_ch1 = brain ph',
                       '!Sample_characteristics_ch1 = brain region',
                       "!Sample_characteristics_ch1 = Cause of death",
                       "!Sample_characteristics_ch1 = disease status",
                       "!Sample_characteristics_ch1 = individual id",
                       "!Sample_characteristics_ch1 = post-mortem interval (in hours)",
                       "!Sample_characteristics_ch1 = rin",
                       "!Sample_characteristics_ch1 = Sex",
                       '!Sample_title',
                       '!Sample_platform_id')]
softData$GSM = rownames(softData)

colnames(softData) = c('deatAge',
                       'ancestry',
                       'brainBank',
                       'pH',
                       'brainRegion',
                       'deathCause',
                       'disease',
                       'invidID',
                       'postMortemInt',
                       'rnaInt',
                       'sex',
                       'title',
                       'platform',
                       'GSM')


# download cel files
dir.create('data-raw/cel/GPL5175')
sapply(softData$GSM,function(x){
    xNow<<-x
    gsmDown(gsm=x,outfile=paste0('data-raw//cel/GPL5175/',x,'.cel'),unzip=F)
})

# attach scandates
softData$scanDate = sapply(softData$GSM, function(x){
    celfileDate(paste0('data-raw/cel/GPL5175/',x, '.cel.gz'))
})
write.design(softData,file='data-raw/TrabzuniRegions/GSE60862_meta.tsv')

 softFile = read.design('data-raw/TrabzuniRegions/GSE60862_meta.tsv')
 softFile = softFile[ !is.na(softFile$pH) 
                      & !softFile$deathCause=='Cancer',]

readOligoCel(softFile$GSM,
             gemmaAnnot= 'data-raw/GemmaAnnots/GPL5175',
             file = 'data-raw/TrabzuniRegions/GSE60862_expression.csv',
             celdir='data-raw/cel//GPL5175')

# outlier detection for trabzuni dataset --------------

humanExp = read.exp('data-raw/TrabzuniRegions/GSE60862_expression.csv')
humanExp = humanExp[!is.na(humanExp$Gene.Symbol),]
humanExp = humanExp[humanExp$Gene.Symbol!='',]
humanExp %<>% mostVariable(threshold=0)

list[genes,expr] = sepExpr(humanExp)
colnames(expr) = gsub('\\.cel\\.gz','',colnames(expr))
softFile = read.design('data-raw/TrabzuniRegions/GSE60862_meta.tsv')

softFile %<>% filter(GSM %in% colnames(expr))

groups = softFile$brainRegion
expression = expr

outliers = sampleMixup(expression, groups)

write.table(outliers, 'data-raw/TrabzuniRegions/outliers', row.names=F,col.names=F,quote=F)

# process trabzuni data (batch correction, removal of outliers) ----------------
humanExp = read.exp('data-raw/TrabzuniRegions//GSE60862_expression.csv')
humanExp = humanExp[!is.na(humanExp$Gene.Symbol),]
humanExp = humanExp[humanExp$Gene.Symbol!='',]
# check expression in all brain regions individually

list[geneData,exprData] = sepExpr(humanExp)
colnames(exprData) = gsub('\\.cel\\.gz','',colnames(exprData))


softFile = read.design('data-raw/TrabzuniRegions/GSE60862_meta.tsv')
softFile = softFile[match(colnames(exprData) ,softFile$GSM),]

medExp = exprData %>%
    unlist %>% median
rowMedians = exprData %>% 
    apply(1, function(x){
        regions =  softFile$brainRegion %>% unique
        sapply(regions, function(y){
            x[softFile$brainRegion %in% y] %>% median
        })
    }) %>% apply(2, max)

humanExp = 
    cbind(geneData, exprData)[rowMedians > medExp,] %>% 
    mostVariable(threshold=0)
list[geneData,exprData] = sepExpr(humanExp)

# remove outliers detected in the previous step
outliers = unlist(read.table('data-raw/TrabzuniRegions/outliers'))
softFile = softFile[!softFile$GSM %in% outliers,]
exprData = exprData[,!colnames(exprData) %in% outliers,]

exprData = as.data.frame(exprData)
rownames(exprData) = geneData$Gene.Symbol

#exprData = exprData[,str_extract(softFile$scanDate,'....-..-..') %in% names(table(str_extract(softFile$scanDate,'....-..-..')))[table(str_extract(softFile$scanDate,'....-..-..'))>1]]
#softFile = softFile[str_extract(softFile$scanDate,'....-..-..') %in% names(table(str_extract(softFile$scanDate,'....-..-..')))[table(str_extract(softFile$scanDate,'....-..-..'))>1],]
set.seed(1)
exprData = ComBat(exprData,batch = kmeans(as.Date(str_extract(softFile$scanDate,'....-..-..')),centers=4)$cluster, mod = model.matrix(~brainRegion,softFile))

#rowMean = exprData %>% apply(1,mean)
#medExp = median(unlist(exprData))
#exprData = exprData[rowMean>medExp,]

write.csv(exprData, file = 'data-raw/TrabzuniRegions/GSE60862_expression_postProcess.csv',quote=F)
write.design(softFile, file = 'data-raw/TrabzuniRegions/GSE60862_meta_postProcess.tsv')

trabzuniRegionsExp = exprData
trabzuniRegionsMeta = softFile
devtools::use_data(trabzuniRegionsExp,
                   trabzuniRegionsMeta, overwrite=TRUE)

# stanley data is taken from Lilah -----------------
load('data-raw/StanleyData/StanleyData.RData')
ogbox::getGemmaAnnot('GPL11532','data-raw/GemmaAnnots/GPL11532',annotType='noParents')

stanleyOut = function(stanleyStud){
    stud = stanleyStud$aned_good
    meta = stanleyStud$Metadata
    
    list[stanleyGenes,stanleyStud] = stud %>% sepExpr
    medExp = median(unlist(stanleyStud))
    rowMedian = stanleyStud %>% apply(1,median)
    
    list[stanleyGenes,stanleyStud] = cbind(stanleyGenes,stanleyStud) %>% 
        mostVariable(.,genes='GeneSymbol',threshold=0) %>% sepExpr
    rownames(stanleyStud) = stanleyGenes$GeneSymbol
    return(list(stanleyStud,meta))   
}

list[stanleyStud1,stanleyMeta1] = stanleyOut(studyFinal$study1)
list[stanleyStud3,stanleyMeta3] = stanleyOut(studyFinal$study3)
list[stanleyStud5,stanleyMeta5] = stanleyOut(studyFinal$study5)
list[stanleyStud7,stanleyMeta7] = stanleyOut(studyFinal$study7)

devtools::use_data(stanleyStud1,
                   stanleyMeta1,
                   stanleyStud3,
                   stanleyMeta3,
                   stanleyStud5,
                   stanleyMeta5,
                   stanleyStud7,
                   stanleyMeta7,
                   overwrite=TRUE)

# Chen (2016) (Sibille) dataset -----------------
ogbox::getGemmaAnnot('GPL11532','data-raw/GemmaAnnots/GPL11532',annotType='noParents')
dir.create('data-raw/ChenCortex/', showWarnings=FALSE)
softDown('GSE71620','data-raw/ChenCortex//GSE71620_family.soft.gz')
system('gunzip data-raw/ChenCortex//GSE71620_family.soft.gz')
softData = softParser('data-raw/ChenCortex/GSE71620_family.soft',expression=FALSE)

#list[softData,expression]  = softParser('data-raw/ChenCortex/GSE71620_family.soft',expression=TRUE)


# expTable = sapply(expression, function(x){
#     x$VALUE
# })

softData = softData[,c('!Sample_characteristics_ch1 = age',
                       '!Sample_characteristics_ch1 = ph',
                       "!Sample_characteristics_ch1 = pmi",
                       '!Sample_characteristics_ch1 = race',
                       "!Sample_characteristics_ch1 = rin",
                       "!Sample_characteristics_ch1 = Sex",
                       "!Sample_characteristics_ch1 = tissue",
                       "!Sample_characteristics_ch1 = tod",
                       "!Sample_title",
                       "!Sample_geo_accession"
)]

names(softData) = c('Age',
                    'pH',
                    'pmi',
                    'race',
                    'rin',
                    'sex',
                    'tissue',
                    'tod',
                    'title',
                    'GSM')
write.design(softData,file='data-raw/ChenCortex/GSE71620_design.tsv')

dir.create('data-raw/cel/GPL11532/',showWarnings=FALSE)

softData = read.design('data-raw/ChenCortex/GSE71620_design.tsv', comment.char='')

sapply(softData$GSM,function(x){
    xNow<<-x
    gsmDown(gsm=x,outfile=paste0('data-raw//cel/GPL11532/',x,'.cel'),unzip=F)
})

expTable = readOligoCel(softData$GSM,
             gemmaAnnot= 'data-raw/GemmaAnnots/GPL11532',
             file = NULL,
             celdir='data-raw/cel/GPL11532')


expTable %<>% filter(Gene.Symbol != '')
expTable %<>% filter(!is.na(Gene.Symbol))


expTable %<>% mostVariable(threshold=expTable %>% sepExpr %>% .[[2]] %>% unlist %>% median,
                           threshFun=median)


write.csv(expTable, 'data-raw/ChenCortex/GSE71620_exp.csv', row.names = F)
chenExpr = read.exp('data-raw/ChenCortex/GSE71620_exp.csv',sep=',')
chenMeta = read.design('data-raw/ChenCortex/GSE71620_design.tsv', comment.char='')

list[chenGenes,chenExpr] = chenExpr %>% sepExpr
rownames(chenExpr) = chenGenes$Gene.Symbol %>% unlist

devtools::use_data(chenExpr,
                  chenMeta,
                  overwrite=TRUE)

# dataset lists ------
brainWholeTissueDatasets = c('trabzuniRegionsExp','chenExpr','stanleyStud1','stanleyStud3','stanleyStud5','stanleyStud7')
brainWholeTissueMetadata = c('trabzuniRegionsMeta','chenMeta', 'stanleyMeta1', 'stanleyMeta3', 'stanleyMeta5','stanleyMeta7')
devtools::use_data(brainWholeTissueDatasets,
                   brainWholeTissueMetadata, overwrite=TRUE)
