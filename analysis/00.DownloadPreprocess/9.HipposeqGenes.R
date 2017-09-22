library(ogbox)
library(devtools)
library(data.table)
library(dplyr)
library(magrittr)
load_all()

dir.create('data-raw/hipposeq',showWarnings = FALSE)
download.file(url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_genes.fpkm_tracking.gz',
              destfile = 'data-raw/hipposeq/hippoSeqExpr.txt.gz')

# download.file(url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_genes.read_group_tracking.txt.gz',
#               destfile = 'data-raw/hipposeq/hippoSeqExpr.txt.gz')

system('gunzip -f data-raw/hipposeq/hippoSeqExpr.txt.gz')

hipposeqExpr = fread('data-raw/hipposeq/hippoSeqExpr.txt',data.table=FALSE)

hipposeqExpr %<>% select(gene_short_name,dg_d_FPKM,dg_v_FPKM,ca4_FPKM,ca3_d_FPKM,ca3_v_FPKM,ca2_FPKM,ca1_d_FPKM,ca1_v_FPKM)

list[gene,expr] = hipposeqExpr %>% sepExpr
expr %<>% {log2(x = .+1)} %>% mutate(dg = (dg_d_FPKM+dg_v_FPKM)/2) %>% select(-dg_d_FPKM,-dg_v_FPKM)


whereAreMyGenes = n_expressoExpr$Gene.Symbol[!n_expressoExpr$Gene.Symbol %in% gene$gene_short_name] %>% mouseSyno(cores = 15)
replaceGenes = whereAreMyGenes %>% sapply(function(x){
    out = gene$gene_short_name[gene$gene_short_name %in% unlist(x)]
    if(length(out)>1){
        out = character(0)
    }
    return(out)
}) %>% unlist
replaceGenes = replaceGenes[!replaceGenes %in% replaceGenes[replaceGenes %>% duplicated]]
replace = names(replaceGenes)
names(replace) = replaceGenes
gene$gene_short_name =gene$gene_short_name %>% replaceElement(replace) %$% newVector




granuleAllowedGenes =gene$gene_short_name[names(expr)[expr %>% apply(1,which.max)] =='dg']


use_data(granuleAllowedGenes,overwrite = TRUE)
