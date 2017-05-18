# comparing dopaminergic cell proportions in substantia nigra samples
library(ogbox)
library(limma)
library(dplyr)
library(magrittr)
library(cowplot)
library(VennDiagram)
library(forcats)
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
# MoranDes %<>% mutate(patient = paste0(Disease, str_extract(Title,pattern='[0-9]+')))
list[geneDat, expDatMoran] = sepExpr(MoranParkinsonsExp)
rownames(expDatMoran) = geneDat$Gene.Symbol

# remove samples from superior frontal gyrus
expDatMoran = expDatMoran[(!MoranDes$Region %in% 'Superior frontal gyrus')]
MoranDes = MoranDes[(!MoranDes$Region %in% 'Superior frontal gyrus') ,]


moranLateralEstimations =  cellTypeEstimate(exprData=cbind(geneDat, expDatMoran[MoranDes$Region %in% 'Lateral substantia nigra']),
                                            genes=genes,
                                            geneColName='Gene.Symbol',
                                            outlierSampleRemove=F,
                                            groups=MoranDes$Disease[MoranDes$Region %in% 'Lateral substantia nigra'],
                                            removeNegatives = T,
                                            PC = 1)

moranMedialEstimation = cellTypeEstimate(exprData=cbind(geneDat, expDatMoran[MoranDes$Region %in% "Medial substantia nigra"]),
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



# Zhang GSE20295 --------
ZhangDes = ZhangParkinsonsMeta %>% filter(brainRegion == 'Whole substantia nigra from postmortem brain')
list[geneDat, expDatZhang] = sepExpr(ZhangParkinsonsExp)
expDatZhang = expDatZhang[ZhangParkinsonsMeta$GSM %in% ZhangDes$GSM]

zhangEstimations =  cellTypeEstimate(exprData=cbind(geneDat, expDatZhang),
                                            genes=genes,
                                            geneColName='Gene.Symbol',
                                            outlierSampleRemove=F,
                                            groups=ZhangDes$diseaseState,
                                            removeNegatives = T,
                                            PC = 1)
ZhangVals = zhangEstimations$estimates %>% sapply(function(x){
    p = wilcox.test(x[ZhangDes$diseaseState == 'Parkinsons disease'], x[ZhangDes$diseaseState == 'Control'])$p.value
})


# dopaminergic estimations plot ------------
geneCounts = c(LesnickEstimations$rotations$Dopaminergic %>% nrow,
               moranLateralEstimations$rotations$Dopaminergic %>% nrow,
               moranMedialEstimation$rotations$Dopaminergic %>% nrow)

LesnickName = paste0('Lesnick\n(n genes = ',LesnickEstimations$rotations$Dopaminergic %>% nrow,')')
MoranLateralName = paste0('Moran Lateral\n(n genes = ',moranLateralEstimations$rotations$Dopaminergic %>% nrow,')')
MoranMediallName = paste0('Moran Medial\n(n genes = ',moranMedialEstimation$rotations$Dopaminergic %>% nrow,')')
ZhangName = paste0('Zhang\n(n genes = ',zhangEstimations$rotations$Dopaminergic %>% nrow,')')

frame1 = data.frame(parkinsons = setNames(c('PD','control'), c(T,F))[LesnickParkinsonsMeta$parkinson %>% as.character],
                    sex = setNames(c('F','M'), c(T,F))[LesnickParkinsonsMeta$female %>% as.character],
                    estimate = scale01(LesnickEstimations$estimates$Dopaminergic),
                    name = LesnickName)

frame2 = data.frame(parkinsons = MoranDes[MoranDes$Region %in% 'Lateral substantia nigra',]$Disease,
                    sex = MoranDes[MoranDes$Region %in% 'Lateral substantia nigra',]$Sex,
                    estimate = scale01(moranLateralEstimations$estimates$Dopaminergic),
                    name = MoranLateralName)

frame3 = data.frame(parkinsons = MoranDes[MoranDes$Region %in% 'Medial substantia nigra',]$Disease,
                    sex = MoranDes[MoranDes$Region %in% 'Medial substantia nigra',]$Sex,
                    estimate = scale01(moranMedialEstimation$estimates$Dopaminergic),
                    name = MoranMediallName)

frame4 = data.frame(parkinsons = ZhangDes$diseaseState %>% replaceElement(c(Control = 'control','Parkinsons disease' = 'PD')) %$% newVector,
                    sex = ZhangDes$gender %>% replaceElement(c('male' ='M',female = 'F')) %$% newVector,
                    estimate = scale01(zhangEstimations$estimates$Dopaminergic),
                    name = ZhangName)

pVals = c(LesnickpVals['Dopaminergic'],
          MoranpValsLateral['Dopaminergic'],
          MoranpValsMedial['Dopaminergic'],
          ZhangVals['Dopaminergic']) %>% signifMarker
signifFrame = data.frame(markers = pVals,
                         x = 1.5,
                         y = 1.0,
                         name = c(LesnickName,
                                  MoranLateralName,
                                  MoranMediallName,
                                  ZhangName))

masterFrame = rbind(frame1,frame2,frame3, frame4)

pEstimate = masterFrame %>%  ggplot(aes( y = estimate, x = parkinsons)) + 
    #geom_point(position= 'jitter',size=3) +
    facet_grid(~name)  +
    theme_cowplot(17) + 
    geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    geom_boxplot(width=0.2,fill = 'lightblue') + 
    # geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 13)) +
    coord_cartesian(ylim = c(-0.10, 1.10)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10)+
    xlab('') +
    ylab('Dopaminergic MGP estimation')

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

dim(LesnickDif)

mm = model.matrix(~ estimate + parkinson,data.frame(parkinson = LesnickParkinsonsMeta$parkinson,
                                                    estimate = LesnickEstimations$estimates$Dopaminergic))
fit <- lmFit(expDatLesnick, mm)
fit <- eBayes(fit)
LesnickDifCorr = topTable(fit, coef=colnames(fit$design)[3],p.value= 0.05,
                    number = Inf)
dim(LesnickDifCorr)

venn = venn.diagram(x=list(Corrected = rn(LesnickDifCorr), Uncorrected = rn(LesnickDif)), filename=NULL)
plot.new()
grid.draw(venn)
# Moran differentially expressed-------
MoranDes = MoranParkinsonsMeta
MoranDes %<>% mutate(patient = paste0(Disease, str_extract(Title,pattern='[0-9]+')))
list[geneDatMoran, expDatMoran] = sepExpr(MoranParkinsonsExp)
rownames(expDatMoran) = geneDatMoran$Gene.Symbol

# remove samples from superior frontal gyrus
expDatMoran = expDatMoran[(!MoranDes$Region %in% 'Superior frontal gyrus')]
MoranDes = MoranDes[(!MoranDes$Region %in% 'Superior frontal gyrus') ,]

mm = model.matrix(~ Disease+ Region,MoranDes)
fit <- lmFit(expDatMoran, mm)
fit <- eBayes(fit)
MoranDif = topTable(fit, coef=colnames(fit$design)[2],
                      #lfc = log(1,base=2),
                      number = Inf,
                      p.value = 0.05
)
dim(MoranDif)

moranEstimations = cellTypeEstimate(exprData=cbind(geneDatMoran, expDatMoran),
                                    genes=genes,
                                    geneColName='Gene.Symbol',
                                    outlierSampleRemove=F,
                                    groups=MoranDes$Disease,
                                    removeNegatives = T,
                                    PC = 1)


mm = model.matrix(~ estimate + parkinson + region,data.frame(parkinson = MoranDes$Disease,
                                                    region = MoranDes$Region,
                                                    estimate = moranEstimations$estimates$Dopaminergic))
fit <- lmFit(expDatMoran, mm)
fit <- eBayes(fit)
MoranDifCorr = topTable(fit, coef=colnames(fit$design)[3],p.value= 0.05,
                          number = Inf)
dim(MoranDifCorr)

venn = venn.diagram(x=list(Corrected = rn(MoranDifCorr), Uncorrected = rn(MoranDif)), filename=NULL)
plot.new()
grid.draw(venn)
# correlations to estimations ----------
# moran PC and estimation
moranPaperGene = data.frame(disease = MoranDes$Disease,
                            paperGenePC = 
                                t(expDatMoran[rownames(expDatMoran) %in% paperGenes,]) %>% 
                                prcomp(.scale=TRUE) %$% x[,1],
                            estimation = moranEstimations$estimates$Dopaminergic)

# overall moran
cor(moranPaperGene$paperGenePC,
    moranPaperGene$estimation, 
    method = 'spearman')
# control moran
moranControlCor = cor(moranPaperGene$paperGenePC[moranPaperGene$disease=='control'],
                      moranPaperGene$estimation[moranPaperGene$disease=='control'],
                      method='spearman') %>% abs
# PD moran
moranPDCor = cor(moranPaperGene$paperGenePC[moranPaperGene$disease=='PD'],
                 moranPaperGene$estimation[moranPaperGene$disease=='PD'],
                 method='spearman') %>% abs


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
lesnickControlCor = cor(lesnickPaperGene$paperGenePC[lesnickPaperGene$disease==FALSE],
                        lesnickPaperGene$estimation[lesnickPaperGene$disease==FALSE],
                        method='spearman') %>% abs
# PD lesnick
lesnickPDCor = cor(lesnickPaperGene$paperGenePC[lesnickPaperGene$disease==TRUE],
                   lesnickPaperGene$estimation[lesnickPaperGene$disease==TRUE],
                   method='spearman') %>% abs

# Zhang paper gene ------------------
zhangPaperGene = data.frame(disease = ZhangDes$diseaseState,
                              paperGenePC = 
                                  t(expDatLesnick[rownames(expDatLesnick) %in% paperGenes,]) %>% 
                                  prcomp(.scale=TRUE) %$% x[,1],
                              estimation = LesnickEstimations$estimates$Dopaminergic)
# overall Zhang
cor(lesnickPaperGene$paperGenePC,
    lesnickPaperGene$estimation, 
    method = 'spearman')
# control Zhang
lesnickControlCor = cor(lesnickPaperGene$paperGenePC[lesnickPaperGene$disease==FALSE],
                        lesnickPaperGene$estimation[lesnickPaperGene$disease==FALSE],
                        method='spearman') %>% abs
# PD Zhang
lesnickPDCor = cor(lesnickPaperGene$paperGenePC[lesnickPaperGene$disease==TRUE],
                   lesnickPaperGene$estimation[lesnickPaperGene$disease==TRUE],
                   method='spearman') %>% abs

# plot's of correlations to gene PC --------------

lesnickPaperGene$disease %<>% ogbox::replaceElement(c('TRUE'='PD','FALSE' = 'control')) %$% newVector
moranPaperGene$paperGenePC = -moranPaperGene$paperGenePC 
masterFrame = rbind(
    data.frame(moranPaperGene,source = 'Moran'),
    data.frame(lesnickPaperGene, source = 'Lesnick'))

annotation = data.frame(source = c('Moran',
                                   'Lesnick'),
                        label = c(paste0('Spearman\'s ρ:\n',
                                       'controls: ',format(lesnickControlCor,digits = 3),'\n',
                                       'PD: ', format(moranPDCor,digits=3)),
                                  paste0('Spearman\'s ρ:\n',
                                         'controls: ',format(moranControlCor,digits = 3),'\n',
                                         'PD: ', format(lesnickPDCor,digits=3))))

genePCestimationPlot = masterFrame %>% ggplot(aes(x =estimation , y = paperGenePC, color = disease)) + 
    facet_grid(~source) +
    theme_cowplot(17) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 17))+
    geom_point(size=3) + 
    stat_smooth(method=lm, se=FALSE) + 
    scale_color_viridis(discrete=TRUE) +
    xlab('Dopaminergic MGP estimation') +
    ylab('PC1 of PD signature gene expression') + 
    geom_text(data = annotation,
              aes(label = label,
                  y = min(masterFrame$paperGenePC),
                  x = max(masterFrame$estimation)),
              color = 'black',
              vjust= 0,
              hjust= 1,
              size = 5)

ggsave(filename = 'analysis//04.MarkerGeneProfiles/publishPlot/genePCestimation.png',
       plot= genePCestimationPlot,
       width=7,height=4.3,units='in')

# plot of correlations to genes -----
# moran
moranGeneCorsPD = cor(moranEstimations$estimates$Dopaminergic[moranPaperGene$disease  == "PD"],
                    expDatMoran[rownames(expDatMoran) %in% paperGenes,moranPaperGene$disease =="PD"] %>% t, 
                    method='spearman') %>% t
moranGeneCorsPD %<>% data.frame(Correlation = ., disease = 'PD', source  ='Moran', gene = rn(.))

moranGeneCorsCont = cor(moranEstimations$estimates$Dopaminergic[moranPaperGene$disease  == 'control'],
                    expDatMoran[rownames(expDatMoran) %in% paperGenes,moranPaperGene$disease  == 'control'] %>% t, 
                    method='spearman') %>% t
moranGeneCorsCont %<>% data.frame(Correlation = ., disease = 'control', source  ='Moran', gene = rn(.))

# lesnick
lesnickGeneCorsCont = cor(LesnickEstimations$estimates$Dopaminergic[lesnickPaperGene$disease  == 'control'],
                        expDatLesnick[rownames(expDatLesnick) %in% paperGenes,lesnickPaperGene$disease  == 'control'] %>% t, 
                        method='spearman') %>% t
lesnickGeneCorsCont %<>% data.frame(Correlation = ., disease = 'control', source  ='Lesnick', gene = rn(.))

lesnickGeneCorsPD = cor(LesnickEstimations$estimates$Dopaminergic[lesnickPaperGene$disease  == 'PD'],
                        expDatLesnick[rownames(expDatLesnick) %in% paperGenes,lesnickPaperGene$disease  == 'PD'] %>% t, 
                          method='spearman') %>% t
lesnickGeneCorsPD %<>% data.frame(Correlation = ., disease = 'PD', source  ='Lesnick', gene = rn(.))

toPlot = rbind(moranGeneCorsCont,
               moranGeneCorsPD,
               lesnickGeneCorsCont,
               lesnickGeneCorsPD)

toPlot$disease %<>% factor(levels=c('control','PD'))

toPlot$gene %<>% fct_reorder(toPlot$Correlation,.desc=TRUE)
    
geneAllestimation = toPlot %>% ggplot(aes(x = gene, y = Correlation, color = disease)) + 
    facet_grid(~source) +
    geom_abline(slope = 0, intercept = 0 ,linetype =2 ) + 
    theme_cowplot(17) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 17))+
    geom_point(size=3, alpha = 0.8) + 
    scale_color_viridis(discrete=TRUE) +
    xlab('') +
    ylab('Dopaminergic MGP-\nGene expression correlation') + 
    theme(axis.text.x  = element_text(angle= 90,vjust = 0.5, size = 13))

ggsave(filename = 'analysis//04.MarkerGeneProfiles/publishPlot/geneAllestimation.png',
       plot= geneAllestimation,
       width=9.5,height=4.3,units='in')

# experimentation on apo-profile relationship --------------
apoSet = list(apo = read.table('geneset.txt',skip = 2)$V1)

LesnickApo = 
    cellTypeEstimate(exprData=LesnickParkinsonsExp,
                     genes=apoSet,
                     geneTransform = NULL,
                     geneColName='Gene.Symbol',
                     outlierSampleRemove=F,
                     groups=setNames(c('parkinson\'s','control'), c(T,F))[LesnickParkinsonsMeta$parkinson %>% as.character],
                     removeNegatives=TRUE,
                     PC = 1)

LesnickFrame = data.frame(source = 'Lesnick',
                          apo = LesnickApo$estimates$apo,
                          dopa = LesnickEstimations$estimates$Dopaminergic,
                          group = replaceElement(LesnickParkinsonsMeta$parkinson, c('TRUE' = 'PD', "FALSE" = 'control'))$newVector)

MoranApo = cellTypeEstimate(exprData=cbind(geneDatMoran, expDatMoran),
                            genes=apoSet,
                            geneTransform = NULL,
                            geneColName='Gene.Symbol',
                            outlierSampleRemove=F,
                            groups=MoranDes$Disease,
                            removeNegatives = T,
                            PC = 1)
moranLateralApo =  cellTypeEstimate(exprData=cbind(geneDat, expDatMoran[MoranDes$Region %in% 'Lateral substantia nigra']),
                                    genes=apoSet,
                                    geneTransform = NULL,
                                    geneColName='Gene.Symbol',
                                    outlierSampleRemove=F,
                                    groups=MoranDes$Disease[MoranDes$Region %in% 'Lateral substantia nigra'],
                                    removeNegatives = T,
                                    PC = 1)
MoranLatFrame  = data.frame(source = 'Moran Lateral',
                            apo = moranLateralApo$estimates$apo,
                            dopa = moranLateralEstimations$estimates$Dopaminergic,
                            group  =MoranDes$Disease[MoranDes$Region %in% 'Lateral substantia nigra'])

moranMedialApo = cellTypeEstimate(exprData=cbind(geneDat, expDatMoran[MoranDes$Region %in% "Medial substantia nigra"]),
                                  genes=apoSet,
                                  geneTransform = NULL,
                                  geneColName='Gene.Symbol',
                                  outlierSampleRemove=F,
                                  groups=MoranDes$Disease[MoranDes$Region %in% "Medial substantia nigra"],
                                  removeNegatives = T,
                                  PC = 1)
MoranMedFrame  = data.frame(source = 'Moran Medial',
                            apo = moranMedialApo$estimates$apo,
                            dopa = moranMedialEstimation$estimates$Dopaminergic,
                            group=MoranDes$Disease[MoranDes$Region %in% "Medial substantia nigra"])


cor(MoranLatFrame$apo[MoranLatFrame$group %in% 'control'],MoranLatFrame$dopa[MoranLatFrame$group %in% 'control'],method = 'spearman')
cor(MoranLatFrame$apo,MoranLatFrame$dopa,method = 'spearman')

rbind(MoranLatFrame,
      MoranMedFrame,
      LesnickFrame) -> hede

hede %>% group_by(source,group) %>% do({cor(.$apo,.$dopa) %>% as.data.frame}) -> yo
names(yo)[3] ='lele'

rbind(MoranLatFrame,
      MoranMedFrame,
      LesnickFrame) %>% ggplot(aes(y = dopa,  x = apo)) + geom_point() + facet_grid(source~group) + geom_smooth(method = 'lm') +
    geom_text(data= yo,x = -9,y = 5,aes(label = format(lele,digits=3)))

