load_all()
library(ggplot2)
library(magrittr)
library(dplyr)

list[geneDat, exp] = n_expressoExpr  %>% 
    filter(!grepl(pattern='\\|',Gene.Symbol)) %>%
    sepExpr
rownames(exp) = geneDat$Gene.Symbol

design = n_expressoSamples

design %<>% filter(PyramidalDeep %in% c('ForebrainCholin',
                                        'ThalamusCholin')) %>%
    mutate(PyramidalDeep = PyramidalDeep %>% factor(levels=c('ForebrainCholin',
                                                             'ThalamusCholin'))) %>% 
    arrange(PyramidalDeep)


cholinergic = exp[geneDat$Gene.Symbol %in% c('Gad1','Gad2','Slc32a1','Slc17a6','Slc18a3'),design$sampleName] %>% 
    as.matrix %>% 
    melt %>% 
    mutate(Var2 = design$ShinyNames[match(Var2,design$sampleName)],
           Var1 = factor(Var1,levels = c('Slc18a3','Gad1','Gad2','Slc32a1','Slc17a6')))

(cholinergic %>% ggplot(aes(x = Var1, y= value, group = Var2, color=Var2))  +
     geom_point(position=position_dodge(width = 0.4), size = 4,fill ='black') +
     scale_x_discrete(name='') + 
     theme_cowplot(17) + 
     theme(legend.position="bottom",
           legend.direction = 'vertical') + 
     geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype='dotted') + 
     scale_y_continuous(name = bquote(log[2]~' expression')) +
     scale_color_manual(name='Cell Type',values = c('darkorange','darkorange4'))) %>%
    ggsave(plot = ., filename='analysis/07.GeneExpressionPlots/cholinergic.png', width=8,height = 4)

# gabaergic
design = n_expressoSamples

design %<>% filter(ShinyNames %in% c('FS Basket (G42)',
                                     'Martinotti (GIN)',
                                     'VIPReln (G30)')) %>%
    mutate(ShinyNames = ShinyNames %>% factor(levels=c('FS Basket (G42)',
                                                       'Martinotti (GIN)',
                                                       'VIPReln (G30)'))) %>% 
    arrange(ShinyNames)  

gabaergic = exp[geneDat$Gene.Symbol %in% c('Sst'),design$sampleName] %>% 
    as.matrix %>% 
    melt %>% 
    mutate(Var2 = design$ShinyNames[match(Var2,design$sampleName)])

(gabaergic %>% ggplot(aes(x= Var2, y = value, color = Var2)) +
     geom_point(size = 4) + theme_cowplot(17) +
     theme(legend.position="none",
           axis.text.x = element_text(angle = 90, hjust = 1)) +
     scale_color_manual(name='Cell Type',values = c('firebrick2','firebrick3','firebrick4')) +
     xlab('') + ylab(bquote('Sst '~log[2]~' expression'))
) %>%  ggsave(plot = ., filename='analysis/07.GeneExpressionPlots/gabaergic.png', width=3,height = 5)

# phani vs chung dopaminergic------------------
design = n_expressoSamples
design %<>% filter(Reference %in% c('Chung et al., 2005',
                                    'Phani et al. 2015')) %>%
    mutate(PyramidalDeep = PyramidalDeep %>% factor(levels=c('Phani et al. 2015',
                                                         'Chung et al., 2005'))) %>% 
    arrange(Reference)


dopaminergic = exp[geneDat$Gene.Symbol %in% c('Map2','Plcb4','Card10','Kifc2'),design$sampleName] %>% 
    as.matrix %>% 
    melt %>% 
    mutate(Var2 = design$Reference[match(Var2,design$sampleName)],
           Var1 = factor(Var1,levels = c('Map2','Plcb4','Card10','Kifc2')))

(dopaminergic %>% ggplot(aes(x = Var1, y= value, group = Var2, color=Var2))  +
     geom_point(position=position_dodge(width = 0.4), size = 4,fill ='black') +
     scale_x_discrete(name='') + 
     theme_cowplot(17) + 
     theme(legend.position="bottom",
           legend.direction = 'vertical') + 
     geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype='dotted') + 
     scale_y_continuous(name = bquote(log[2]~' expression')) +
     scale_color_manual(name='Cell Type',values = c('black','gray'), guide = guide_legend(title = ""))) %>%
    ggsave(plot = ., filename='analysis/07.GeneExpressionPlots/dopaminergic.png', width=5,height = 4)

# Olig1 and Fam114a1 in cortex----------------
design = n_expressoSamples

regionSamples = memoReg(design = design,regionNames = 'Region',groupNames = 'ShinyNames',regionHierarchy=regionHierarchy)

design %<>% filter(!is.na(regionSamples$Cortex_ShinyNames)) %>% 
    arrange(MajorType,Neurotransmitter,PyramidalDeep)

cortical = exp[geneDat$Gene.Symbol %in% c('Olig1','Fam114a1'),design$sampleName] %>% as.matrix %>% melt %>% 
    mutate(Var2 = design$ShinyNames[match(Var2,design$sampleName)] %>% factor(levels=design$ShinyNames %>% unique))

(cortical %>% ggplot(aes(x = Var2, y = value, color = Var2)) + 
     geom_point(size = 1) +
     theme_cowplot(8) + 
     xlab('') +
     ylab(bquote(log[2]~' expression'))+
     theme(axis.text.x = element_text(angle = 90, hjust = 1),
           panel.background = element_rect(fill=NA, color="black",size=0.3,linetype='solid'))+
     facet_grid(Var1~.,scales='free',switch='y') + 
     scale_color_manual(values = cellColors(),guide = FALSE)) %>%
    ggsave(plot = .,filename = 'analysis/07.GeneExpressionPlots/olig1_Fam114a1.png', width = 5.25,height = 9.5,units='cm')

# Ddc in oligodendrocyte --------------
design = n_expressoSamples
regionSamples = memoReg(design = design,regionNames = 'Region',groupNames = 'ShinyNames',regionHierarchy=regionHierarchy)
design %<>% filter((!is.na(regionSamples$Cortex_ShinyNames))| design$ShinyNames %in% 'Dopaminergic') %>% 
    arrange(MajorType,Neurotransmitter,PyramidalDeep)

cortical = exp[geneDat$Gene.Symbol %in% c('Ddc'),design$sampleName] %>% as.matrix %>% melt %>% 
    mutate(Source = design$Reference[match(Var2,design$sampleName)],
           Var2 = design$ShinyNames[match(Var2,design$sampleName)] %>% factor(levels=design$ShinyNames %>% unique))           

cortical$Groups = cortical %>% apply(1,function(x){
    if(x['Var2'] == 'Oligodendrocyte'){
        return(x['Source'])
    } else if(x['Var2'] == 'Dopaminergic'){
        return('dopaminergic')
    } else{
        return('others')
    }
})

cortical$Groups %<>% factor(levels=c('Cahoy et al., 2008', 'Doyle et al., 2008', 'Fomchenko et al 2011', 'dopaminergic', 'others'))

(cortical %>% ggplot(aes(x = Groups, y = value, color = Groups)) + 
     geom_point(size = 4) +
     theme_cowplot(17) + 
     xlab('') +
     ylab(bquote('Ddc '~log[2]~' expression'))+
     theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
     scale_color_viridis(discrete =  TRUE, guide=FALSE)) %>% ggsave(plot = .,filename = 'analysis/07.GeneExpressionPlots/Ddc.png', width = 3.5,height = 6)


