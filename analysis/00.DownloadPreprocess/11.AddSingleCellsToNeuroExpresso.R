library(ogbox)
load_all()

singleCells = ogbox::read.design('data-raw/Mouse_Cell_Type_Data/singleCellMatchings.tsv')
tasicUnassigned =  (TasicMouseMeta$primary_type %>% unique)[! (TasicMouseMeta$primary_type %>% unique) %in% (singleCells$Tasic %>% str_split(', ') %>% unlist %>% trimElement(''))]

TasicPrimaryMean = TasicMouseMeta$primary_type %>% unique %>% lapply(function(x){
    TasicMouseExp[,TasicMouseMeta$sample_title[TasicMouseMeta$primary_type %in% x]] %>% apply(1,mean)
}) %>% as.data.frame
names(TasicPrimaryMean)  =  TasicMouseMeta$primary_type %>% unique


# filter expression values (low level filter)
TasicPrimaryMean = TasicPrimaryMean[(TasicPrimaryMean %>% apply(1,max))>(TasicPrimaryMean %>% unlist %>% median),]



TasicPrimaryMeanComparable = TasicPrimaryMean %>% 
    apply(2,qNormToValues,values =  n_expressoExpr %>%
              sepExpr %>% {.[[2]]} %>% unlist) %>%
    as.df
rownames(TasicPrimaryMeanComparable) = rn(TasicPrimaryMean)



TasicPrimaryMeanSubset =  TasicPrimaryMean %>% {.[rn(.) %in% n_expressoExpr$Gene.Symbol,]}
n_ExpressoSubset = n_expressoExpr[match(rn(TasicPrimaryMeanSubset), n_expressoExpr$Gene.Symbol),]

TasicPrimaryMeanComparableRows = mapply(function(tasic,neuro){
    qNormToValues(tasic,neuro)
    },
    t(TasicPrimaryMeanSubset) %>% as.df ,n_ExpressoSubset %>% sepExpr %>% {.[[2]]} %>% t %>% as.df) %>% t %>% as.df

names(TasicPrimaryMeanComparableRows) = names(TasicPrimaryMean)



sampleLines = singleCells$Tasic %>% str_split(', ')

samples = unlist(sampleLines)


meltedSingleCells = samples %>% sapply(function(x){
    if(!x==''){
        findInList(x , sampleLines)
    } else{
        NULL
    }
}) %>% unlist %>% {
    out = singleCells[.,]
    out$sample = names(.)
    return(out)
}


set.seed(1)
meltedSingleCells %<>%
    {data.frame(sampleName =.$sample,
               originalIndex = stri_rand_strings(nrow(.),6),
               GSE =  'GSE71585',
               samples = .$Tasic,
               MajorType = .$MajorType,
               Neurotransmitter = .$Neurotransmitter,
               ShinyNames = .$ShinyNames  %>% replaceElement(NA,'') %$%newVector,
               PyramidalDeep = .$PyramidalDeep %>% replaceElement(NA,'') %$%newVector,
               BroadTypes = NA,
               DopaSelect = NA,
               Description = NA,
               Age = NA,
               Region = .$Region,
               Platform = 'RNAseq',
               Method = 'RNA-seq',
               Reference = 'Tasic et al.',
               PMID = 26727548,
               RegionToChildren = TRUE,
               RegionToParent = TRUE,
               Normalize = TRUE,
               Normalize2= TRUE,
               Notes = '',
               stringsAsFactors = FALSE)}

TasicPrimaryMean = TasicPrimaryMean[meltedSingleCells$sampleName]
TasicPrimaryMeanComparable = TasicPrimaryMeanComparable[meltedSingleCells$sampleName]

use_data(meltedSingleCells, overwrite = TRUE)
use_data(TasicPrimaryMean, overwrite = TRUE)
use_data(TasicPrimaryMeanComparable, overwrite = TRUE)
write.design(meltedSingleCells,file = 'data-raw/Mouse_Cell_Type_Data/meltedSingleCells.tsv')
write.csv(data.frame(Gene.Symbol = rn(TasicPrimaryMeanComparable),TasicPrimaryMeanComparable,check.names = FALSE),file = 'data-raw/Mouse_Cell_Type_Data/TasicPrimaryMeanComparable.csv',row.names=FALSE )
write.csv(data.frame(Gene.Symbol = rn(TasicPrimaryMean),TasicPrimaryMean,check.names=FALSE),file = 'data-raw/Mouse_Cell_Type_Data/TasicPrimaryMean.csv',row.names=FALSE )


n_expressoSamplesWithRNAseq = rbind(n_expressoSamples, meltedSingleCells)

n_expressoExprWithRNAseq = cbind(n_expressoExpr[n_expressoExpr$Gene.Symbol %in% rn(TasicPrimaryMeanComparable),],
                                 TasicPrimaryMeanComparable[match(n_expressoExpr$Gene.Symbol[n_expressoExpr$Gene.Symbol %in% rn(TasicPrimaryMeanComparable)],
                                                                  rn(TasicPrimaryMeanComparable)),])

write.design(n_expressoSamplesWithRNAseq,'data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq.tsv')
use_data(n_expressoSamplesWithRNAseq,overwrite = TRUE)
use_data(n_expressoExprWithRNAseq,overwrite = TRUE)
write.csv(n_expressoExprWithRNAseq,file = 'data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq.csv',row.names=FALSE )


n_expressoExprWithRNAseqRowNorm = cbind(n_expressoExpr[n_expressoExpr$Gene.Symbol %in% rn(TasicPrimaryMeanComparableRows),],
                                        TasicPrimaryMeanComparableRows[match(n_expressoExpr$Gene.Symbol[n_expressoExpr$Gene.Symbol %in% rn(TasicPrimaryMeanComparableRows)],
                                                                         rn(TasicPrimaryMeanComparableRows)),])

write.csv(n_expressoExprWithRNAseqRowNorm,file = 'data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseqRowNorm.csv',row.names=FALSE )
use_data(n_expressoExprWithRNAseqRowNorm,overwrite = TRUE)


n_expressoSamples2 = read.design('data-raw/Mouse_Cell_Type_Data/n_expressoSamples2.tsv')
n_expressoExpr2 = read.exp("data-raw/Mouse_Cell_Type_Data/n_expressoExpr2.csv")
n_expressoSamplesWithRNAseq2 = rbind(n_expressoSamples2, meltedSingleCells)
write.design(n_expressoSamplesWithRNAseq2,'data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq2.tsv')


n_expressoExprWithRNAseq2 = cbind(n_expressoExpr2[n_expressoExpr2$Gene.Symbol %in% rn(TasicPrimaryMeanComparable),],
                                 TasicPrimaryMeanComparable[match(n_expressoExpr2$Gene.Symbol[n_expressoExpr2$Gene.Symbol %in% rn(TasicPrimaryMeanComparable)],
                                                                  rn(TasicPrimaryMeanComparable)),])
write.csv(n_expressoExprWithRNAseq2,file = 'data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq2.csv',row.names=FALSE )


n_expressoExprWithRNAseq2RowNorm = cbind(n_expressoExpr2[n_expressoExpr2$Gene.Symbol %in% rn(TasicPrimaryMeanComparableRows),],
                                         TasicPrimaryMeanComparableRows[match(n_expressoExpr2$Gene.Symbol[n_expressoExpr2$Gene.Symbol %in% rn(TasicPrimaryMeanComparableRows)],
                                                                   rn(TasicPrimaryMeanComparableRows)),])

write.csv(n_expressoExprWithRNAseq2RowNorm,file = 'data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq2RowNorm.csv',row.names=FALSE )
