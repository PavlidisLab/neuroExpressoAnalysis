library(limma)
cels = paste0('data-raw//cel//GPL1261/', c('GSM1621133','GSM1621134','GSM1621135',
                                           'GSM1621130','GSM1621131','GSM1621132'),'.cel')




affy = ReadAffy(filenames = cels)
affy = affy::rma(affy)

affy = gemmaAnnot(affy, paste0('data-raw/GemmaAnnots//','/','GPL1261'))


affy %<>% mostVariable(threshold=affy %>% sepExpr %>% .[[2]] %>% unlist %>% median,
                       threshFun=median)

list[genes,exp] = affy %>% sepExpr
rownames(exp) = genes$Gene.Symbol
astrocytesReactive = exp
colnames(astrocytesReactive) = c('astro1',
                                 'astro2',
                                 'astro3',
                                 'reactive1',
                                 'reactive2',
                                 'reactive3')

dir.create('data-raw/reactiveAstros',showWarnings=FALSE)
write.csv(astrocytesReactive,'data-raw/reactiveAstros/humanRNASeq.csv',quote=F)
devtools::use_data(astrocytesReactive,overwrite=TRUE)   
