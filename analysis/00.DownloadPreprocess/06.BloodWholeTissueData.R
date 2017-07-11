library(XLConnect)
devtools::load_all()

downloadData = FALSE
downloadGemma = TRUE

#lymphoma dataset ----------
dir.create('data-raw/wholeBloodDatasets',showWarnings=FALSE)
if(downloadGemma){
    getGemmaAnnot('GPL570','data-raw/GemmaAnnots/GPL570',annotType='noParents',overwrite=TRUE)
}
if(downloadData){
    gseDown(GSE='GSE65135',regex='lymph',outDir='data-raw/cel/GPL570')
}
cels =celFiles('data-raw/cel/GPL570',full.names=T) 
cels = cels[grep(regexMerge(gsmFind('GSE65135', 'lymph')), cels)]
affy = ReadAffy(filenames = cels)
norm = affy::rma(affy)
annotated = gemmaAnnot(norm, 'data-raw/GemmaAnnots/GPL570')
names(annotated) = gsub('[.]cel','',names(annotated))
annotated = mostVariable(allDataPre=annotated,
                         threshold=annotated %>% sepExpr %>% .[[2]] %>% unlist %>% median,
                         threshFun=median)

LYMPHexpr = annotated

write.csv(annotated, paste0('data-raw/wholeBloodDatasets/','lymphoma.csv'), row.names = FALSE)

devtools::use_data(LYMPHexpr, overwrite=TRUE)

# pbmc dataset ---------------
data = httr::GET('https://cibersort.stanford.edu/inc/inc.download.page.handler.php?file=PBMCs-Fig3a-HumanHT-12-V4.txt') %>% httr::content('text')
PBMCs = read.table(textConnection(data),sep='\t',header=T)
write.table(PBMCs, 'data-raw/wholeBloodDatasets/PBMCs.csv', sep=',', quote=FALSE, row.names=F)
PBMCexpr = PBMCs
devtools::use_data(PBMCexpr, overwrite=TRUE)



# cell counts ----
download.file('http://www.nature.com/nmeth/journal/v12/n5/source_data/nmeth.3337-f3.xlsx',
              destfile='data-raw/wholeBloodDatasets/fig3.xlsx')
cibersort = loadWorkbook('data-raw/wholeBloodDatasets/fig3.xlsx')
PBMCcibersort = readWorksheet(cibersort, sheet = 1, header = TRUE)
LYMPHCibersort = readWorksheet(cibersort, sheet = 3, header = TRUE)

gsms = sapply(LYMPHCibersort$Sample.ID, function(x){gsmFind(GSE='GSE65135', regex=x)})
LYMPHCibersort$Sample.ID = gsms
names(LYMPHCibersort) = c('Sample.ID', 'Cell.type', 'RealCounts', 'CibersortPred')

names(PBMCcibersort) = c('Sample.ID', 'Cell.type', 'RealCounts', 'CibersortPred')


write.table(LYMPHCibersort, 'data-raw/wholeBloodDatasets/LYMPHCibersort',sep = '\t', quote=FALSE, row.names=FALSE)
write.table(PBMCcibersort, 'data-raw/wholeBloodDatasets/PBMCcibersort.tsv',sep = '\t', quote=FALSE, row.names=FALSE)
devtools::use_data(PBMCcibersort,
                   LYMPHCibersort, overwrite=TRUE)

