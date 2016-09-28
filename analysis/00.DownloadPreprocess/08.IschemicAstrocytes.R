library(ogbox)
library(magrittr)
library(affy)
library(stringr)
devtools::load_all()


gseDown(GSE='GSE35338',
        outDir='data-raw//cel//GPL1261')



cels = paste0('data-raw//cel//GPL1261/',
              ogbox::gsmFind(GSE='GSE35338'),
                             '.cel')

# library(limma)
# cels = paste0('data-raw//cel//GPL1261/', c('GSM1621133','GSM1621134','GSM1621135',
#                                            'GSM1621130','GSM1621131','GSM1621132'),'.cel')
# 
# 
affy = ReadAffy(filenames = cels)
affy = affy::rma(affy)

affy = gemmaAnnot(affy, paste0('data-raw/GemmaAnnots//','/','GPL1261'))


affy %<>% mostVariable(threshold=affy %>% sepExpr %>% .[[2]] %>% unlist %>% median,
                       threshFun=median)

list[genes,exp] = affy %>% sepExpr
rownames(exp) = genes$Gene.Symbol
astrocytesIschemic = exp


colnames(astrocytesIschemic) %<>% gsub(pattern='[.]cel',replacement='',.) %>% sapply(ogbox::geoTitle)

    
dir.create('data-raw/ischemicAstros',showWarnings=FALSE)
write.csv(astrocytesIschemic,'data-raw/ischemicAstros/ischemicAstros.csv',quote=F)
devtools::use_data(astrocytesIschemic,overwrite=TRUE)   


# get the differentially expressed genes between groups -------
groups = cn(astrocytesIschemic) %>% str_extract(pattern=regexMerge(c('MCAO','LPS', 'saline','sham')))

# no lps as we don't really care about it here
comparisons = combn(c('MCAO', 'saline','sham'),2)

difGenes = lapply(1:ncol(comparisons), function(i){
    model = groups[groups %in% comparisons[,i]]
    data = astrocytesIschemic[,groups %in% comparisons[,i]]
    mm = model.matrix(~ model,as.df(model))
    fit <- lmFit(data, mm)
    fit <- eBayes(fit)
    difGenes = topTable(fit, coef=colnames(fit$design)[2],lfc=log(2,base=2),
                          #lfc = log(1,base=2),
                          number = Inf, 
                          p.value = 0.05)
    rn(difGenes)    
})

names(difGenes) = sapply(1:ncol(comparisons), function(i){
    paste(comparisons[,i],collapse=' ')
})


library(VennDiagram)
venn = venn.diagram(difGenes,filename=NULL)
plot.new()
grid.draw(venn)


# look at the genes that are differentially expressed between sham surgery and MCAO
ischemiaGenes = c(difGenes$`MCAO sham`)
devtools::use_data(ischemiaGenes,overwrite=TRUE)


