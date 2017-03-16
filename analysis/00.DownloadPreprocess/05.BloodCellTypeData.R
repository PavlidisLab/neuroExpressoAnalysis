# download blood cell type datasets
library(XLConnect)
devtools::load_all()

downloadData = FALSE
downloadGemma = TRUE

# human blood------
dir.create('data-raw/GemmaAnnots',showWarnings = FALSE)
if(downloadGemma){
    ogbox::getGemmaAnnot(chipName='GPL96', chipFile='data-raw/GemmaAnnots/GPL96',annotType='noParents',overwrite = TRUE)
}

dir.create('data-raw/HumanBloodCellTypeData', showWarnings=FALSE)

if(downloadData){
    download.file('http://www.nature.com/nmeth/journal/v12/n5/extref/nmeth.3337-S2.xls',
                  destfile='data-raw/HumanBloodCellTypeData/supp2.xls')
}
bloodDes = loadWorkbook('data-raw/HumanBloodCellTypeData/supp2.xls')
bloodDes = readWorksheet(bloodDes, sheet = 2, header = TRUE)

# special case because reasons....
bloodDes$Sample.ID[bloodDes$Sample.ID %in% 'A_MF_2hrEosinophils_U133A'] = 'A_MF_2hrEosinophils'

bloodDes$lm11 = c(rep('B Cells',15),rep('PCs', 7), rep('CD8 T cells', 4), rep('CD4', 14), rep('GammaDeltaT', 2),
                    rep('NK', 15), rep('MonoMacro', 30), rep('Dendritic', 12), rep('Mast', 4), rep('Eos', 2), rep('Neutrophils',8))

bloodDes$originalIndex = as.numeric(factor(bloodDes$Abreviated.name))
names(bloodDes)[names(bloodDes) %in% 'Sample.ID'] ='sampleName'
names(bloodDes)[names(bloodDes) %in% 'Abreviated.name'] ='lm22'

write.design(bloodDes, 'data-raw/HumanBloodCellTypeData/humanBloodCellsSamples.tsv')
humanBloodCellsSamples = bloodDes
devtools::use_data(humanBloodCellsSamples,overwrite=TRUE)

gsms = bloodDes$sampleName[grepl('GSM',bloodDes$sampleName)]

if(downloadData){
    sapply(gsms,function(gsm){
        gsmDown(gsm, paste0('data-raw/cel/GPL96/',gsm,'.CEL'))
    })
    
    # download the ones that not GSMs
    links = c('http://www.ebi.ac.uk/arrayexpress/files/E-MEXP-750/E-MEXP-750.raw.1.zip/TN_U133A_1.CEL',
              'http://www.ebi.ac.uk/arrayexpress/files/E-MEXP-750/E-MEXP-750.raw.1.zip/TN_U133A_2.CEL',
              'http://www.ebi.ac.uk/arrayexpress/files/E-MEXP-750/E-MEXP-750.raw.1.zip/TN_U133A_3.CEL',
              'http://www.ebi.ac.uk/arrayexpress/files/E-MEXP-750/E-MEXP-750.raw.1.zip/CXCR5hiICOShi_U133A_1.CEL',
              'http://www.ebi.ac.uk/arrayexpress/files/E-MEXP-750/E-MEXP-750.raw.1.zip/CXCR5hiICOShi_U133A_2.CEL',
              'http://www.ebi.ac.uk/arrayexpress/files/E-MEXP-750/E-MEXP-750.raw.1.zip/CXCR5hiICOShi_U133A_3.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/Tcells/gamma%20delta/A_TS_RN_gdTcells_U133A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/Tcells/gamma%20delta/A_TS_RN_gdTcellsREP_A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/mastcell/Control/exp2/A_LW_mastcellctrl_U133A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/mastcell/Control/exp1/A_MF_ControlMASTCELL_U133A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/mastcell/IgE/exp1/A_MF_IgEMASTCELL_U133A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/mastcell/IgE/exp2/A_LW_mastcellIgE_U133A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/eosinophils/pma/A_MF_2hrEosinophils.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/eosinophils/control/A_MF_ControlEosinophil.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/neutrophils/control/exp2/A_LW_neutrophil_U133A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/neutrophils/control/exp1/A_MF_neutrophils_U133A.CEL',
              'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/neutrophils/LPS/A_TS_MSNeutroLPS_U133A.CEL')
    
    links %>% sapply(function(x){
        if(file.exists(paste0('data-raw//cel//GPL96/',basename(x)))){
            print('You already have this file. Skipping')
            return(NULL)
        }
        
        download.file(x, paste0('data-raw//cel//GPL96/',basename(x)))
    })
}
cels = paste0('data-raw/cel/GPL96/',bloodDes$sampleName,'.CEL')
affy = ReadAffy(filenames = cels)
norm = affy::rma(affy)
annotated = gemmaAnnot(norm, 'data-raw/GemmaAnnots/GPL96')
names(annotated) = gsub('[.]CEL','',names(annotated))
#write.csv(annotated, 'data-raw/HumanBloodCellTypeData/bloodCellsExp.csv', row.names = F)
annotated = quantileNorm(annotated)

#humanBloodCellsExpALLPROBES = annotated
#devtools::use_data(humanBloodCellsExpALLPROBES,overwrite=TRUE)

dir.create('data-raw/HumanBloodCellTypeData/',showWarnings=FALSE)

humanBloodCellsExp= mostVariableCT(annotated,
                              'data-raw/HumanBloodCellTypeData/humanBloodCellsExp.csv',
                              cellTypeColumn = 'lm22',
                              design=bloodDes)
devtools::use_data(humanBloodCellsExp,overwrite=TRUE)

# mouse blood ------
if(downloadGemma){
    ogbox::getGemmaAnnot(chipName='GPL1261', chipFile='data-raw/GemmaAnnots/GPL96',annotType='noParents',overwrite = TRUE)
}
mouseBloodDes = read.design('data-raw/MouseBloodCellTypeData/mouseBloodDes.tsv')
# here for historic reasons. delete eventually
# mouseHumanDictionary22 = c(Treg = 'Tregs',
#                          Neutrophil = "PMNs",
#                          Monoctye = 'Monos',
#                          'T cells follicular helper' = 'Tfh cells',
#                          "NK activated" = 'NK cells+',
#                          "NK resting" = "NK cells-",
#                          "B cell naive" = "Naïve B cells",
#                          "Macrophage M0" = "M0-MΦs",
#                          "Macrophage M1" ="M1-MΦs" ,
#                          "Macrophage M2" = "M2-MΦs",
#                          "T_CD4_naive" = "CD4 naïve T cells",
#                          "T_CD8" = "CD8 T cells",
#                          "T_CD4_memory" = "CD4 memory T cells-",
#                          Eosinophil = 'Eos',
#                          "DC resting" ="DCs-" ,
#                          "DC activated" =  "DCs+",
#                          "T_CD4_activated" = "CD4 memory T cells+",
#                          "B cell memory" = "Memory B cells",
#                          "T cell gamma delta" = "γδ T cells",
#                          "Plasma cell" = "PCs",
#                          "Mast cell resting" = "MCs-",
#                          "Mast cell activated" = "MCs+")
# 
# mouseHumanDictionary11 = c(CD4 = "CD4",
#                            Neutrophil = 'Neutrophils',
#                            Macrophage = 'MonoMacro',
#                            NK = 'NK',
#                            "B Cell" = "B Cells",
#                            CD8 = "CD8 T cells",
#                            Eosinophil = 'Eos',
#                            Dendritic = 'Dendritic',
#                            "T cell gamma delta"=  'GammaDeltaT',
#                            'Plasma cell' = 'PCs',
#                            "Mast cell" = 'Mast')
# mouseBloodDes$lm22 %<>% replaceElement(mouseHumanDictionary22) %$%newVector
# mouseBloodDes$lm11 %<>% replaceElement(mouseHumanDictionary11) %$% newVector
# write.design(mouseBloodDes, file='data-raw/MouseBloodCellTypeData/mouseBloodDes.tsv')
gsms = mouseBloodDes$GSM
if(downloadData){
    sapply(gsms,function(gsm){
        gsmDown(gsm, paste0('data-raw/cel/GPL1261/',gsm,'.cel'))
    })
}
cels = paste0('data-raw/cel/GPL1261/',gsms,'.cel')
affy = ReadAffy(filenames = cels)
norm = affy::rma(affy)
annotated = gemmaAnnot(norm, 'data-raw/GemmaAnnots/GPL1261')
names(annotated) = gsub('[.]cel','',names(annotated))
#write.csv(annotated, 'data-raw/HumanBloodCellTypeData/bloodCellsExp.csv', row.names = F)
annotated = quantileNorm(annotated)

#mouseBloodCellsExpALLPROBES = annotated 
#devtools::use_data(mouseBloodCellsExpALLPROBES,overwrite=TRUE)


mouseBloodCellsExp = mostVariableCT(annotated,
                                   'data-raw/MouseBloodCellTypeData/mouseBloodCellsExp.csv',
                                   cellTypeColumn = 'lm22',
                                   design=mouseBloodDes)
mouseBloodCellsSamples = mouseBloodDes
devtools::use_data(mouseBloodCellsSamples,overwrite=TRUE)
devtools::use_data(mouseBloodCellsExp,overwrite=TRUE)

