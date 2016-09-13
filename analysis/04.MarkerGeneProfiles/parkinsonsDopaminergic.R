# comparing dopaminergic cell proportions in substantia nigra samples
library(ogbox)
library(limma)
library(dplyr)
library(magrittr)
devtools::load_all()


# Lesnick GSE7621--------

genes = mouseMarkerGenes$Midbrain[
    !grepl('Microglia_',names(mouseMarkerGenes$Midbrain))]


LesnickEstimations = 
    cellTypeEstimate(exprData=LesnickParkinsonsExp,
                     genes=genes,
                     geneColName='Gene.Symbol',
                     outlierSampleRemove=F,
                     groups=setNames(c('parkinson\'s','control'), c(T,F))[LesnickParkinsonsMeta$parkinson %>% as.character],
                     removeNegatives=TRUE,
                     PC = 1)

LesnickpVals = LesnickEstimations$estimates %>% sapply(function(x){
    p = wilcox.test(x[LesnickParkinsonsMeta$parkinson == TRUE], x[LesnickParkinsonsMeta$parkinson == FALSE])$p.value
})

# Moran GSE8397 ----------
MoranDes = MoranParkinsonsMeta
MoranDes %<>% mutate(patient = paste0(Disease, str_extract(Title,pattern='[0-9]+')))
list[geneDat, expDat] = sepExpr(MoranParkinsonsExp)
rownames(expDat) = geneDat$Gene.Symbol

# remove samples from superior frontal gyrus
expDat = expDat[(!MoranDes$Region %in% 'Superior frontal gyrus')]
MoranDes = MoranDes[(!MoranDes$Region %in% 'Superior frontal gyrus') ,]


moranLateralEstimations =  cellTypeEstimate(exprData=cbind(geneDat, expDat[MoranDes$Region %in% 'Lateral substantia nigra']),
                                            genes=genes,
                                            geneColName='Gene.Symbol',
                                            outlierSampleRemove=F,
                                            groups=MoranDes$Disease[MoranDes$Region %in% 'Lateral substantia nigra'],
                                            removeNegatives = T,
                                            PC = 1)

moranMedialEstimation = cellTypeEstimate(exprData=cbind(geneDat, expDat[MoranDes$Region %in% "Medial substantia nigra"]),
                                         genes=genes,
                                         geneColName='Gene.Symbol',
                                         outlierSampleRemove=F,
                                         groups=MoranDes$Disease[MoranDes$Region %in% "Medial substantia nigra"],
                                         removeNegatives = T,
                                         PC = 1)


MoranpValsLateral = moranLateralEstimations$estimates %>% sapply(function(x){
    MoranDes = MoranDes[MoranDes$Region %in% 'Lateral substantia nigra',]
    p = wilcox.test(x[MoranDes$Disease == 'PD'], x[MoranDes$Disease == 'control'])$p.value
})

MoranpValsMedial = moranMedialEstimation$estimates %>% sapply(function(x){
    MoranDes = MoranDes[MoranDes$Region %in% 'Medial substantia nigra',]
    p = wilcox.test(x[MoranDes$Disease == 'PD'], x[MoranDes$Disease == 'control'])$p.value
})


# dopaminergic estimations plot ------------
frame1 = data.frame(parkinsons = setNames(c('PD','control'), c(T,F))[LesnickParkinsonsMeta$parkinson %>% as.character],
                    sex = setNames(c('F','M'), c(T,F))[LesnickParkinsonsMeta$female %>% as.character],
                    estimate = scale01(LesnickEstimations$estimates$Dopaminergic),
                    name = 'Lesnick\nGSE7621')

frame2 = data.frame(parkinsons = MoranDes[MoranDes$Region %in% 'Lateral substantia nigra',]$Disease,
                    sex = MoranDes[MoranDes$Region %in% 'Lateral substantia nigra',]$Sex,
                    estimate = scale01(moranLateralEstimations$estimates$Dopaminergic),
                    name = 'Moran\nGSE8397 lateral')

frame3 = data.frame(parkinsons = MoranDes[MoranDes$Region %in% 'Medial substantia nigra',]$Disease,
                    sex = MoranDes[MoranDes$Region %in% 'Medial substantia nigra',]$Sex,
                    estimate = scale01(moranMedialEstimation$estimates$Dopaminergic),
                    name = 'Moran\nGSE8397 medial')

pVals = c(LesnickpVals['Dopaminergic'],
          MoranpValsLateral['Dopaminergic'],
          MoranpValsMedial['Dopaminergic']) %>% signifMarker
signifFrame = data.frame(markers = pVals,
                         x = 1.5,
                         y = 1.0,
                         name = c('Lesnick\nGSE7621',
                                  'Moran\nGSE8397 lateral',
                                  'Moran\nGSE8397 medial'))

masterFrame = rbind(frame1,frame2,frame3)

pEstimate = masterFrame %>%  ggplot(aes( y = estimate, x = parkinsons)) + 
    #geom_point(position= 'jitter',size=3) +
    facet_grid(~name)  +
    theme_cowplot(17) + 
    geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    geom_boxplot(width=0.2,fill = 'lightblue') + 
    # geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 14)) +
    coord_cartesian(ylim = c(-0.10, 1.10)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10)+
    xlab('') +
    ylab('Cell type profile estimation')

ggsave(plot=pEstimate,
       filename='analysis//04.MarkerGeneProfiles/publishPlot/dopaminergicEstimation.png',width=7,height=4.5,units='in')

# correlations to differentially expressed genes ---------------

# gene list published in the paper
c('CDC42', 'FGF13',
  "HSPB1","SNCA",
  "MKNK2", "TF",
  "AMPH", "BEX1",
  "JMJD6", "NSF",
  "SUB1", "SV2B",
  "SYT1", "SNAP25",
  "STMN2", "RGS4",
  "SNX10", "PRKAR2B",
  "NEFL", "MDH1",
  "CHGB",
  "NFASC") %>% sort -> paperGenes


# Lesnick differentially expressed -----
list[geneDatLesnick,expDatLesnick] = sepExpr(LesnickParkinsonsExp)
rownames(expDatLesnick) = geneDatLesnick$Gene.Symbol

mm = model.matrix(~ parkinson,LesnickParkinsonsMeta)
fit <- lmFit(expDatLesnick, mm)
fit <- eBayes(fit)
LesnickDif = topTable(fit, coef=colnames(fit$design)[2],
               #lfc = log(1,base=2),
               number = Inf, 
               p.value = 0.05
               )

# Moran differentially expressed-------
MoranDes = MoranParkinsonsMeta
MoranDes %<>% mutate(patient = paste0(Disease, str_extract(Title,pattern='[0-9]+')))
list[geneDatMoran, expDatMoran] = sepExpr(MoranParkinsonsExp)
rownames(expDatMoran) = geneDatMoran$Gene.Symbol

# remove samples from superior frontal gyrus
expDatMoran = expDatMoran[(!MoranDes$Region %in% 'Superior frontal gyrus')]
MoranDes = MoranDes[(!MoranDes$Region %in% 'Superior frontal gyrus') ,]

mm = model.matrix(~ Disease,MoranDes)
fit <- lmFit(expDatMoran, mm)
fit <- eBayes(fit)
MoranDif = topTable(fit, coef=colnames(fit$design)[2],
                      #lfc = log(1,base=2),
                      number = Inf, 
                      p.value = 0.05
)


moranEstimations = cellTypeEstimate(exprData=cbind(geneDatMoran, expDatMoran),
                                    genes=genes,
                                    geneColName='Gene.Symbol',
                                    outlierSampleRemove=F,
                                    groups=MoranDes$Disease,
                                    removeNegatives = T,
                                    PC = 1)


# correlations to estimations ----------
# moran PC and estimation
moranPaperGene = data.frame(disease = MoranDes$Disease,
                            paperGenePC = 
                                t(expDatMoran[rownames(expDatMoran) %in% paperGenes,]) %>% 
                                prcomp(.scale=TRUE) %$% x[,1],
                            estimation = moranEstimations$estimates$Dopaminergic)
# plot moran
moranPaperGene %>% ggplot(aes(
    ))

# overall moran
cor(moranPaperGene$paperGenePC,
    moranPaperGene$estimation, 
    method = 'spearman')
# control moran
cor(moranPaperGene$paperGenePC[moranPaperGene$disease=='control'],
    moranPaperGene$estimation[moranPaperGene$disease=='control'],
    method='spearman')
# PD moran
cor(moranPaperGene$paperGenePC[moranPaperGene$disease=='PD'],
    moranPaperGene$estimation[moranPaperGene$disease=='PD'],
    method='spearman')


# lesnick PC and estimation
lesnickPaperGene = data.frame(disease = LesnickParkinsonsMeta$parkinson,
                            paperGenePC = 
                                t(expDatLesnick[rownames(expDatLesnick) %in% paperGenes,]) %>% 
                                prcomp(.scale=TRUE) %$% x[,1],
                            estimation = LesnickEstimations$estimates$Dopaminergic)
# overall lesnick
cor(lesnickPaperGene$paperGenePC,
    lesnickPaperGene$estimation, 
    method = 'spearman')
# control lesnick
cor(lesnickPaperGene$paperGenePC[lesnickPaperGene$disease==FALSE],
    lesnickPaperGene$estimation[lesnickPaperGene$disease==FALSE],
    method='spearman')
# PD lesnick
cor(lesnickPaperGene$paperGenePC[lesnickPaperGene$disease==TRUE],
    lesnickPaperGene$estimation[lesnickPaperGene$disease==TRUE],
    method='spearman')


# mean corelation of positive genes moran
t(expDatMoran[rn(expDatMoran) %in% paperGenes,
              MoranDes$Disease %in% 'PD']) %>%
    cor(moranEstimations$estimates$Dopaminergic[MoranDes$Disease %in% 'PD'],
        method='spearman') %>% 
    as.vector %>% {.[.>0]} %>%
    mean
t(expDatMoran[rn(expDatMoran) %in% paperGenes,
              MoranDes$Disease %in% 'control']) %>%
    cor(moranEstimations$estimates$Dopaminergic[MoranDes$Disease %in% 'control'], 
        method='spearman') %>% 
    as.vector %>% {.[.>0]} %>%
    mean

# mean corelation of positive genes lesnick
t(expDatLesnick[rn(expDatLesnick) %in% paperGenes,
                LesnickParkinsonsMeta$parkinson %in% TRUE]) %>%
    cor(LesnickEstimations$estimates$Dopaminergic[LesnickParkinsonsMeta$parkinson %in% TRUE],
        method='spearman') %>% 
    as.vector %>% {.[.>0]} %>% mean

t(expDatLesnick[rn(expDatLesnick) %in% paperGenes,
                LesnickParkinsonsMeta$parkinson %in% FALSE]) %>%
    cor(LesnickEstimations$estimates$Dopaminergic[LesnickParkinsonsMeta$parkinson %in% FALSE],
        method='spearman') %>% 
    as.vector %>% {.[.>0]} %>% mean


t(expDatMoran[rownames(expDatMoran) %in% paperGenes,MoranDes$Disease %in% 'PD']) %>%
    cor(moranEstimations$estimates$Dopaminergic[MoranDes$Disease %in% 'PD'], method='spearman') %>% 
    as.vector %>% {.[.<0]} %>% mean
t(expDatMoran[rownames(expDatMoran) %in% paperGenes,MoranDes$Disease %in% 'control']) %>%
    cor(moranEstimations$estimates$Dopaminergic[MoranDes$Disease %in% 'control'], method='spearman') %>% 
    as.vector %>% {.[.<0]} %>% mean


# expression of positive genes in dopaminergic samples
posGenes = t(expDatMoran[rn(expDatMoran) %in% paperGenes,MoranDes$Disease %in% 'control']) %>%
    cor(moranEstimations$estimates$Dopaminergic[MoranDes$Disease %in% 'control'], method='spearman') %>% 
    as.vector %>% 
{rn(expDatMoran)[rn(expDatMoran) %in% paperGenes][. > 0]}

negGenes =  t(expDatMoran[rn(expDatMoran) %in% paperGenes,MoranDes$Disease %in% 'control']) %>%
    cor(moranEstimations$estimates$Dopaminergic[MoranDes$Disease %in% 'control'], method='spearman') %>% 
    as.vector %>% 
{rn(expDatMoran)[rn(expDatMoran) %in% paperGenes][. < 0]}

list[n_exprGenes,n_expr] = n_expressoExpr %>% sepExpr

n_expr[n_exprGenes$Gene.Symbol %in% (posGenes %>% human2mouse %$% mouseGene),
       n_expressoStudies$PyramidalDeep %in% 'Dopaminergic'] %>% apply(1,mean)

n_expr[n_exprGenes$Gene.Symbol %in% (negGenes %>% human2mouse %$% mouseGene),
       n_expressoStudies$PyramidalDeep %in% 'Dopaminergic'] %>% apply(1,mean)




t(expDatMoran[rownames(expDatMoran) %in% rn(MoranDif),]) %>% prcomp(.scale=TRUE) %$% x[,1] %>% 
    plot(moranEstimations$estimates$Dopaminergic,col= toColor(MoranDes$Disease,c('PD' = 'red', 'control' = 'black'))$cols)

cor(t(expDatMoran[rownames(expDatMoran) %in% rn(MoranDif),]) %>% prcomp(.scale=TRUE) %$% x[,1],
    moranEstimations$estimates$Dopaminergic, method = 'spearman')

t(expDatMoran[rownames(expDatMoran) %in% rn(MoranDif),]) %>%
    cor(moranEstimations$estimates$Dopaminergic) %>% as.vector %>% plot



t(expDatLesnick[rownames(expDatLesnick) %in% paperGenes,]) %>% prcomp(.scale=TRUE) %$% x[,1] %>% 
    plot(LesnickEstimations$estimates$Dopaminergic,col= toColor(MoranDes$Disease,c('PD' = 'red', 'control' = 'black'))$cols)

cor(t(expDatLesnick[rownames(expDatLesnick) %in% paperGenes,]) %>% prcomp(.scale=TRUE) %$% x[,1],
    LesnickEstimations$estimates$Dopaminergic, method = 'spearman')

t(expDatLesnick[rownames(expDatLesnick) %in% paperGenes,]) %>%
    cor(LesnickEstimations$estimates$Dopaminergic) %>% as.vector %>% plot


t(expDatLesnick[rownames(expDatLesnick) %in% rn(LesnickDif),]) %>% prcomp(.scale=TRUE) %$% x[,1] %>% 
    plot(LesnickEstimations$estimates$Dopaminergic,col= toColor(MoranDes$Disease,c('PD' = 'red', 'control' = 'black'))$cols)

cor(t(expDatLesnick[rownames(expDatLesnick) %in% rn(LesnickDif),]) %>% prcomp(.scale=TRUE) %$% x[,1],
    LesnickEstimations$estimates$Dopaminergic, method = 'spearman')

t(expDatLesnick[rownames(expDatLesnick) %in% rn(LesnickDif),]) %>%
    cor(LesnickEstimations$estimates$Dopaminergic) %>% as.vector %>% plot
