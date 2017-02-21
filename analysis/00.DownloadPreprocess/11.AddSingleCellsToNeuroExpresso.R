singleCells = ogbox::read.design('data-raw/Mouse_Cell_Type_Data/singleCellMatchings.tsv')
tasicUnassigned =  (TasicMouseMeta$primary_type %>% unique)[! (TasicMouseMeta$primary_type %>% unique) %in% (singleCells$Tasic %>% str_split(', ') %>% unlist %>% trimElement(''))]

TasicPrimaryMean = TasicMouseMeta$primary_type %>% unique %>% lapply(function(x){
    TasicMouseExp[,TasicMouseMeta$sample_title[TasicMouseMeta$primary_type %in% x]] %>% apply(1,mean)
}) %>% as.data.frame
names(TasicPrimaryMean)  =  TasicMouseMeta$primary_type %>% unique



TasicPrimaryMeanComparable = TasicPrimaryMean %>%  
    qNormToValues(values = n_expressoExpr %>% unlist)%>%
    as.data.frame
rownames(TasicPrimaryMeanComparable) =TasicPrimaryMean %>% rn 

singleCells = singleCells %>% select(TasicBased,PyramidalDeep,ShinyNames,Tasic)


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
               MajorType = .$TasicBased %>%
                   sapply(function(x){
                       n_expressoSamples$MajorType[which(n_expressoSamples$TasicBased %in% x)[1]]
                       }),
               Neurotransmitter = .$TasicBased %>%
                   sapply(function(x){
                       n_expressoSamples$Neurotransmitter[which(n_expressoSamples$TasicBased %in% x)[1]]
                   }),
               ShinyNames = .$ShinyNames  %>% replaceElement(NA,'') %$%newVector,
               TasicBased = .$TasicBased  %>% replaceElement(NA,'') %$%newVector,
               PyramidalDeep = .$PyramidalDeep %>% replaceElement(NA,'') %$%newVector,
               BroadTypes = NA,
               DopaSelect = NA,
               Description = NA,
               Age = NA,
               Region = 'Cortex',
               Platform = 'TasicRNAseq',
               Method = 'RNA-seq',
               Reference = 'Tasic et al.',
               PMID = 26727548,
               RegionToChildren = TRUE,
               RegionToParent = TRUE,
               Normalize = TRUE,
               Normalize2= TRUE,
               Notes = '',
               stringsAsFactors = FALSE)}

use_data(meltedSingleCells, overwrite = TRUE)
use_data(TasicPrimaryMean, overwrite = TRUE)
use_data(TasicPrimaryMeanComparable, overwrite = TRUE)

n_expressoSamplesWithRNAseq = rbind(n_expressoSamples, meltedSingleCells)

n_expressoExprWithRNAseq = cbind(n_expressoExpr[n_expressoExpr$Gene.Symbol %in% rn(TasicPrimaryMeanComparable),],
                                 TasicPrimaryMeanComparable[match(n_expressoExpr$Gene.Symbol[n_expressoExpr$Gene.Symbol %in% rn(TasicPrimaryMeanComparable)],
                                                                  rn(TasicPrimaryMeanComparable)),])

write.design(n_expressoSamplesWithRNAseq,'data-raw/Mouse_Cell_Type_Data/meltedWithRnaSeq.tsv')
use_data(n_expressoSamplesWithRNAseq,overwrite = TRUE)
use_data(n_expressoExprWithRNAseq,overwrite = TRUE)
write.csv(n_expressoExprWithRNAseq,file = 'data-raw/Mouse_Cell_Type_Data/finalExpWithRnaSeq.tsv',row.names=FALSE )
