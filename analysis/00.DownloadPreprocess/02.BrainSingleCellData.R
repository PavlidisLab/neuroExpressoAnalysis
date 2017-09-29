library(XLConnect)
library(magrittr)
library(ogbox)
library(dplyr)

# mouse rna seq data download Zeisel et al. download------
dir.create('data-raw/ZeiselMouse', showWarnings=FALSE)
download.file(url= 'http://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt', 
              destfile='data-raw/ZeiselMouse/mouseRNASeq_Zeisel 2015.txt')
download.file(url = 'http://science.sciencemag.org/highwire/filestream/628248/field_highwire_adjunct_files/1/aaa1934_TableS1.xlsx',
              destfile = 'data-raw/ZeiselMouse/markerGenes.xlsx')
download.file(url = 'http://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt',
              destfile = 'data-raw/ZeiselMouse/expression_spikes_17-Aug-2014.txt')


# mouse rna seq data download Zeisel et al. process------
rnaSeq = read.table('data-raw/ZeiselMouse/mouseRNASeq_Zeisel 2015.txt', sep= '\t', comment.char= "",stringsAsFactors=F)
rnaMeta = rnaSeq[1:10,3:ncol(rnaSeq)]
rnaMeta = as.data.frame(t(rnaMeta))
colnames(rnaMeta) = rnaSeq[1:10,2]

rnaExp = rnaSeq[12:nrow(rnaSeq),3:ncol(rnaSeq)]
rnaExp = apply(rnaExp,2,as.numeric)
rnaExp = matrix(unlist(rnaExp), nrow = nrow(rnaExp))
rownames(rnaExp) = rnaSeq[12:nrow(rnaSeq),1]
# rnaCelIDs = as.numeric(as.character(rnaSeq[12:nrow(rnaSeq), 2]))
rnaExp = rnaExp[,rnaMeta$tissue %in% 'sscortex']
rnaMeta = rnaMeta[rnaMeta$tissue %in% 'sscortex',]
# remove low expressed ones 
maximExp = apply(rnaExp,1,max)
rnaExp = rnaExp[maximExp >= 1,]

rnaMeta$tissue %<>% as.character
rnaMeta$`group #` %<>% as.numeric
rnaMeta$`total mRNA mol` %<>% as.numeric
rnaMeta$well %<>% as.numeric
rnaMeta$sex %<>% as.numeric
rnaMeta$age %<>% as.numeric
rnaMeta$diameter %<>% as.numeric
rnaMeta$cell_id %<>% as.char
rnaMeta$level1class %<>% as.char
rnaMeta$level2class %<>% as.char

ZeiselMouseExp = rnaExp
ZeiselMouseMeta = rnaMeta


devtools::use_data(ZeiselMouseExp,overwrite=TRUE)
metadata is not needed in this study
devtools::use_data(ZeiselMouseMeta,overwrite=TRUE)

# tasic et al. allen institute single cell data download -------------
dir.create('data-raw/TasicMouse', showWarnings=FALSE)
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71585/suppl/GSE71585_RefSeq_RPKM.csv.gz',
              'data-raw/TasicMouse/GSE71585_RefSeq_RPKM.csv.gz')
system('gunzip data-raw/TasicMouse/GSE71585_RefSeq_RPKM.csv.gz')

download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71585/suppl/GSE71585_Clustering_Results.csv.gz',
              'data-raw/TasicMouse/GSE71585_Clustering_Results.csv.gz')
system('gunzip data-raw/TasicMouse/GSE71585_Clustering_Results.csv.gz')


allenBrain = read.exp('data-raw/TasicMouse/GSE71585_RefSeq_RPKM.csv')
genes = allenBrain$gene
allenBrain = allenBrain[,-1]
rownames(allenBrain) = genes
allenBrain %<>% as.matrix


# fix the order
allenMeta = read.csv('data-raw/TasicMouse/GSE71585_Clustering_Results.csv',stringsAsFactors=FALSE)
allenBrain = allenBrain[,allenMeta$sample_title]

maxExp = allenBrain %>% apply(1,max)
allenBrain = allenBrain[maxExp>0,]


TasicMouseExp = allenBrain
TasicMouseMeta = allenMeta
devtools::use_data(TasicMouseExp,overwrite=TRUE)
devtools::use_data(TasicMouseMeta,overwrite=TRUE)



# human rna seq data download from darmanis download -----
gsms = gsmFind('GSE67835')
dir.create('data-raw/DarmanisHuman/raw',recursive=TRUE showWarnings=FALSE)

sapply(gsms, function(gsm){
    print(gsm)
    page = getURL(paste0('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',gsm))
    fileURL = URLdecode(str_extract(page,'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM.*?csv%2Egz'))
    if (len(fileURL) == 0){
        if (warnings){
            warning(paste(gsm,"doesn't have a file attached"))
        }
        return(invisible(F))
    }
    download.file(fileURL,paste0('data-raw/DarmanisHuman/raw/',gsm,'.csv.gz'))
    system(paste0('gunzip -f "',paste0('data-raw/DarmanisHuman/raw/',gsm,'.csv.gz'),'"'))  
})

files = list.files('data-raw/DarmanisHuman/raw', full.names=T)

allExpr = sapply(files,function(x){
    read.table(x, sep='\t')[,2]
})
print('files read')

singleGenes = read.table(files[1], sep='\t', stringsAsFactors=FALSE)[,1]
singleGenes = sapply(singleGenes, trimWS)
rownames(allExpr) = singleGenes
colnames(allExpr) = gsub('[.]csv','',basename(colnames(allExpr)))
write.csv(allExpr,'data-raw/DarmanisHuman/humanRNASeq.csv',quote=FALSE)


softDown('GSE67835',file='data-raw/DarmanisHuman/humanRNAseq.soft.gz')
system('gunzip data-raw/DarmanisHuman/humanRNAseq.soft.gz')
humanMeta = softParser('data-raw/DarmanisHuman/humanRNAseq.soft',expression=FALSE)
humanMeta = humanMeta[c('!Sample_characteristics_ch1 = cell type',
                        '!Sample_characteristics_ch1 = age',
                        '!Sample_characteristics_ch1 = experiment_sample_name',
                        '!Sample_characteristics_ch1 = c1 chip id')]
names(humanMeta) = c('cellType','age','sample','chip')
humanMeta$GSM = rownames(humanMeta)
write.design(humanMeta,'data-raw/DarmanisHuman/humanRNASeq_metadat.tsv')

# human rna seq data download from darmanis process -----
rnaExp = read.table('data-raw/DarmanisHuman/humanRNASeq.csv', sep= ',', comment.char= "",stringsAsFactors=F, row.names=1, header=T)
rnaExp = rnaExp[1:(nrow(rnaExp)-3),]
maximExp = apply(rnaExp,1,max)
rnaExp = rnaExp[maximExp >= 1,]
genes = rn(rnaExp)
rnaExp = apply(rnaExp,2,as.numeric)
rownames(rnaExp) = genes

rnaMeta = read.design('data-raw/DarmanisHuman/humanRNASeq_metadat.tsv')

DarmanisHumanExp = rnaExp
DarmanisHumanMeta = rnaMeta
devtools::use_data(DarmanisHumanExp,overwrite=TRUE)
devtools::use_data(DarmanisHumanMeta,overwrite=TRUE)


singleCellRNASeqDatasets = c('ZeiselMouseExp', 'TasicMouseExp', 'DarmanisHumanExp')
devtools::use_data(singleCellRNASeqDatasets,overwrite=TRUE)
