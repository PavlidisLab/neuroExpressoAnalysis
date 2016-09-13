# download blood cell type datasets
library(XLConnect)

# human blood------
ogbox::getGemmaAnnot(chipName='GPL96', chipFile='data-raw/GemmaAnnots/GPL96',annotType='noParents')
dir.create('data-raw/HumanBloodCellTypeData', showWarnings=FALSE)

download.file('http://www.nature.com/nmeth/journal/v12/n5/extref/nmeth.3337-S2.xls',
              destfile='data-raw/HumanBloodCellTypeData/supp2.xls')
bloodDes = loadWorkbook('data-raw/HumanBloodCellTypeData/supp2.xls')
bloodDes = readWorksheet(bloodDes, sheet = 2, header = TRUE)

# special case because reasons....
bloodDes$Sample.ID[bloodDes$Sample.ID %in% 'A_MF_2hrEosinophils_U133A'] = 'A_MF_2hrEosinophils'

bloodDes$Leuk11 = c(rep('B Cells',15),rep('PCs', 7), rep('CD8 T cells', 4), rep('CD4', 14), rep('GammaDeltaT', 2),
                    rep('NK', 15), rep('MonoMacro', 30), rep('Dendritic', 12), rep('Mast', 4), rep('Eos', 2), rep('Neutrophils',8))

bloodDes$originalIndex = as.numeric(factor(bloodDes$Abreviated.name))
names(bloodDes)[names(bloodDes) %in% 'Sample.ID'] ='sampleName'
write.design(bloodDes, 'data-raw/HumanBloodCellTypeData/humanBloodCellsSamples.tsv')
humanBloodCellsSamples = bloodDes
devtools::use_data(humanBloodCellsSamples,overwrite=TRUE)

gsms = bloodDes$sampleName[grepl('GSM',bloodDes$sampleName)]

sapply(gsms,function(gsm){
    gsmDown(gsm, paste0('data-raw/cel/GPL96/',gsm,'.CEL'))
})


cels = paste0('data-raw/cel/GPL96/',bloodDes$sampleName,'.CEL')
affy = ReadAffy(filenames = cels)
norm = rma(affy)
annotated = gemmaAnnot(norm, 'data-raw/GemmaAnnots/GPL96')
names(annotated) = gsub('[.]CEL','',names(annotated))
#write.csv(annotated, 'data-raw/HumanBloodCellTypeData/bloodCellsExp.csv', row.names = F)
annotated = quantileNorm(annotated)

humanBloodCellsExp= mostVariableCT(annotated,
                              'data-raw/HumanBloodCellTypeData/humanBloodCellsExp.csv',
                              cellTypeColumn = 'Abreviated.name',
                              design=bloodCellsSamples)
devtools::use_data(humanBloodCellsExp,overwrite=TRUE)

