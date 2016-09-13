library(limma)
library(data.table)
library(stringr)
library(sva)
devtools::load_all()

# trabzuni dataset ----------------------
ogbox::getGemmaAnnot('GPL5175','data-raw/GemmaAnnots/GPL5175',annotType='noParents')

dir.create('data-raw/TrabzuniRegions/',showWarnings=FALSE)
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60862/soft/GSE60862_family.soft.gz',
              destfile='data-raw/TrabzuniRegions/GSE60862_family.soft.gz')
system('gunzip data-raw/TrabzuniRegions/SE60862_family.soft.gz')

# deal with softfile
softData = softParser('data-raw/TrabzuniRegions/GSE60862_family.soft')

list[softData,expression]  = softParser('data-raw/TrabzuniRegions/GSE60862_family.soft',expression=TRUE)

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
write.design(softData,file='data-raw/TrabzuniRegions/GSE60862_meta')
# softData = read.design('data/GSE60862_meta')
dir.create('data/cel/GPL5175')
sapply(softData$GSM,function(x){
    xNow<<-x
    gsmDown(gsm=x,outfile=paste0('data//cel/GPL5175/',x),unzip=F)
})

# attach scandates
softData$scanDate = sapply(softData$GSM, function(x){
    celfileDate(paste0('/home/omancarci/masterOfCellTypes/cel/GPL5175/',x, '.cel.gz'))
})
write.design(softData,file='data-raw/TrabzuniRegions/GSE60862_meta.tsv')

softFile = read.design('data-raw/TrabzuniRegions/GSE60862_meta.tsv')
softFile = softFile[ !is.na(softFile$pH) 
                     & !softFile$deathCause=='Cancer',]

readHumanCel(softFile$GSM,'data-raw/TrabzuniRegions/GSE60862_expression.csv',humanDir='data/cel//GPL5175')

# outlier detection for trabzuni dataset --------------

humanExp = fread('data-raw/TrabzuniRegions/GSE60862_expression.csv')
humanExp = humanExp[!is.na(humanExp$Gene.Symbol),]
humanExp = humanExp[humanExp$Gene.Symbol!='',]

list[genes,expr] = sepExpr(humanExp)
colnames(expr) = gsub('\\.cel\\.gz','',colnames(expr))
softFile = read.design('data-raw/TrabzuniRegions/GSE60862_meta.tsv')

groups = unique(softFile$brainRegion)
pairwise = combn(groups,2)

difs = vector(mode = 'list', length = ncol(pairwise))

# get differentially expressed genes between any two regions
for (i in 1:ncol(pairwise)){
    subsetExpr = expr[,names(expr) %in% softFile$GSM[softFile$brainRegion %in% pairwise[,i]],with=F]
    subsetExpr = as.matrix(subsetExpr)
    rownames(subsetExpr) = genes$Gene.Symbol
    
    subsetGroups = softFile$brainRegion[match(colnames(subsetExpr), softFile$GSM)]  
    
    mm = model.matrix(~ subsetGroups,data.frame(subsetGroups))
    fit <- lmFit(subsetExpr, mm)
    fit <- eBayes(fit)
    dif = topTable(fit, coef=colnames(fit$design)[2],
                   lfc = log(5,base=2),
                   number = Inf, 
                   p.value = 1e-5)
    if (nrow(dif)==0){
        print(paste('no genes found for',paste(pairwise[,i],collapse='/') ,'pair'))
        next
    }
    difs[[i]] = dif$ID
    print(i)
}

# detect outliers based on differentially expressed genes inside a single group
outliers = lapply(1:len(groups), function(i){
    relevant = unlist(difs[apply(pairwise,2,function(x){
        groups[[i]] %in%  x  
    })])
    subsetExpr = expr[,names(expr) %in% softFile$GSM[softFile$brainRegion %in% groups[i]],with=F]
    subsetExpr = subsetExpr[genes$Gene.Symbol %in% relevant,]
    pca = prcomp(t(subsetExpr))
    names(boxplot.stats(pca$x[,1])$out)
})

write.table(unlist(outliers),'data-raw/TrabzuniRegions/outliers', row.names=F,col.names=F,quote=F)

# process trabzuni data (batch correction, removal of outliers) ----------------
humanExp = fread('data-raw/TrabzuniRegions//GSE60862_expression.csv')
humanExp = humanExp[!is.na(humanExp$Gene.Symbol),]
humanExp = humanExp[humanExp$Gene.Symbol!='',]

humanExp = mostVariable(humanExp,treshold=0) 
list[geneData,exprData] = sepExpr(humanExp)
# get rid of file extensions
colnames(exprData) = str_extract(string=colnames(exprData), pattern='GSM.*?(?=\\.)')


softFile = read.design('data-raw/TrabzuniRegions/GSE60862_meta.tsv')
softFile = softFile[match(colnames(exprData) ,softFile$GSM),]

# remove outliers detected in the previous step
outliers = unlist(read.table('data-raw/TrabzuniRegions/outliers'))
softFile = softFile[!softFile$GSM %in% outliers,]
exprData = exprData[,!colnames(exprData) %in% outliers,with=F]

# remove white matter from the mix
# exprData = exprData[,!softFile$brainRegion %in% 'white matter',with=F]
# softFile = softFile[!softFile$brainRegion %in% 'white matter',]

exprData = as.data.frame(exprData)
rownames(exprData) = geneData$Gene.Symbol

#exprData = exprData[,str_extract(softFile$scanDate,'....-..-..') %in% names(table(str_extract(softFile$scanDate,'....-..-..')))[table(str_extract(softFile$scanDate,'....-..-..'))>1]]
#softFile = softFile[str_extract(softFile$scanDate,'....-..-..') %in% names(table(str_extract(softFile$scanDate,'....-..-..')))[table(str_extract(softFile$scanDate,'....-..-..'))>1],]
set.seed(1)
exprData = ComBat(exprData,batch = kmeans(as.Date(str_extract(softFile$scanDate,'....-..-..')),centers=4)$cluster, mod = model.matrix(~brainRegion,softFile))
rowmax = exprData %>% apply(1,max)
medExp = median(unlist(exprData))

exprData = exprData[rowmax>medExp,]

write.csv(exprData, file = 'data-raw/TrabzuniRegions/GSE60862_expression_postProcess.csv',quote=F)
write.design(softFile, file = 'data-raw/TrabzuniRegions/GSE60862_meta_postProcess.tsv')

trabzuniRegionsExp = exprData
trabzuniRegionsMeta = softFile
devtools::use_data(trabzuniRegionsExp,
                   trabzuniRegionsMeta)

# stanley data is taken from Lilah -----------------
load('data-raw/StanleyData/StanleyData.RData')
ogbox::getGemmaAnnot('GPL11532','data-raw/GemmaAnnots/GPL11532',annotType='noParents')

stanleyOut = function(stanleyStud){
    stud = stanleyStud$aned_good
    meta = stanleyStud$Metadata
    list[stanleyGenes,stanleyStud] = stud %>%
        mostVariable(.,genes='GeneSymbol',treshold=.[,-1] %>% unlist %>% median) %>% 
        sepExpr
    
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
dir.create('data-raw/ChenCortex/', showWarnings=FALSE)
softDown('GSE71620','data-raw/ChenCortex//GSE71620_family.soft.gz')
system('gunzip data-raw/ChenCortex//GSE71620_family.soft.gz')
list[softData,expression]  = softParser('data-raw/ChenCortex/GSE71620_family.soft',expression=TRUE)


expTable = sapply(expression, function(x){
    x$VALUE
})

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
write.design(softData,file='data-raw/ChenCortex/GSE71620_design')

tableGenes = gemmaGeneMatch(probesets=expression[[1]]$ID_REF,chipFile='data-raw/GemmaAnnots/GPL11532')

expTable = data.frame(Gene.Symbol = tableGenes, expTable)
expTable %<>% filter(Gene.Symbol != '')
expTable %<>% mostVariable(treshold=expTable[,-1] %>% unlist %>% median)

write.csv(expTable, paste0('data-raw/ChenCortex/','GSE71620_exp.csv'), row.names = F, quote = F)

chenExpr = read.exp('data-raw/ChenCortex/GSE71620_exp.csv')
chenMeta = read.design('data-raw/ChenCortex/GSE71620_design', comment.char='')

list[chenGenes,chenExpr] = chenExpr %>% sepExpr
rownames(chenExpr) = chenGenes %>% unlist

devtools::use_data(chenExpr,
                  chenMeta,
                  overwrite=TRUE)

# dataset lists ------
brainWholeTissueDatasets = c('trabzuniRegionsExp','chenExpr','stanleyStud1','stanleyStud3','stanleyStud5','stanleyStud7')
brainWholeTissueMetadata = c('trabzuniRegionsMeta','chenMeta', 'stanleyMeta1', 'stanleyMeta3', 'stanleyMeta5','stanleyMeta7')
devtools::use_data(brainWholeTissueDatasets,
                   brainWholeTissueMetadata)
