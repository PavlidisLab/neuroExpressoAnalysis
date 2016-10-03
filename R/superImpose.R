#' @export
superImpose = function(genes, expression, CIBERSORTpred, geneSymbol, outDir,PC=1){
    theirPred = CIBERSORTpred
    table = {if(!is.null(outDir)){
        paste0(outDir,'/rotTable_', names(genes))
    }else {
        NULL}
    }
    estimates = cellTypeEstimate(exprData = expression, 
                                 genes= genes,
                                 geneColName = geneSymbol,
                                 outlierSampleRemove = F,
                                 synonymTaxID = NULL,
                                 geneTransform = NULL, 
                                 groups = NA,
                                 tableOut = table,
                                 indivGenePlot = NULL,
                                 seekConsensus = F,
                                 plotType = NULL,
                                 PC = PC)
    
    dir.create(outDir,showWarnings=FALSE)
    correlations = vector(mode='list',length = len(estimates$estimates))
    for (i in 1:len(estimates$estimates)){
        if(is.na(estimates$estimates[[i]][1])){
            next
        }
        print(i)
        counts = theirPred %>%
            filter(Cell.type == names(estimates$estimates[i])) %>%
            select(RealCounts)
        cybersort = theirPred %>%
            filter(Cell.type == names(estimates$estimates[i])) %>%
            select(CibersortPred)
        
        estimates$estimates[[i]] = scale01(estimates$estimates[[i]])
        
        scaleCyber = scaleToInt(cybersort,max(estimates$estimates[[i]]), min(estimates$estimates[[i]]))
        
        frame1 = data.frame(estimates = estimates$estimates[[i]], counts, scaleCyber)
        frame2 = data.frame(cybersort, counts)
        
        
        frame = data.frame(estimates = estimates$estimates[[i]],cybersort ,counts)
        
        # simplify the graph by reducing the number of breaks, x axis of second plot is ignored
        countRange = frame1$RealCounts %>% range
        xbreak1 = round((countRange[2] - countRange[1])/2 + countRange[1],digits=-1)
        xbreak2 = round(countRange[2],digits=-1)
        
        if(xbreak1 < countRange[1]){
            xbreak1 = countRange[1] %>% ceiling
        }
        if(xbreak2>countRange[2]){
            xbreak2 = countRange[2] %>% floor
        }
        
        correlations[[i]] = c(CIBERSORT = format(cor(frame1$CibersortPred,frame1$RealCounts,method='spearman'),digits=3),
                              "Marker Genes" = format(cor(frame1$estimates,frame1$RealCounts,method='spearman'),digits=3))
        
        p1 = ggplot(frame1, aes(y= estimates,x = RealCounts)) + geom_point(size = 6, shape = 15) +  theme_bw() +
            scale_x_continuous(breaks = c(xbreak1,xbreak2)) + 
            scale_y_continuous(breaks = c(0.5, 1),limits=c(0,1.2)) + 
            xlab('') +
            ylab('') +
            geom_segment(aes(x=RealCounts,
                             y = estimates,
                             xend = RealCounts,
                             yend = CibersortPred)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            annotate('text', y= 1.2,x= min(frame1$RealCounts),vjust = 1, hjust = 0, size = 11,
                     label = paste0('spearman correlation:\n',
                                    'Marker genes: ', format(cor(frame1$estimates,frame1$RealCounts,method='spearman'),digits=3), '\n',
                                    'Cibersort: ',format(cor(frame1$CibersortPred,frame1$RealCounts,method='spearman'),digits=3))
            ) + geom_line(size=3,alpha=1/10) +
            ggtitle(paste0(names(estimates$estimates[i]),'\n(n genes = ',
                    estimates$rotations[[i]] %>% nrow, ')')) + 
            theme(axis.text= element_text(size=30),
                  plot.title = element_text(size=35)) 
        
        cibersortRange = frame2$CibersortPred  %>% range
        #ybreak1 = round(cibersortRange[1],digits=-1) + (mod(cibersortRange[1],5) <5)*10
        ybreak1 = round((cibersortRange[2] - cibersortRange[1])/2 + cibersortRange[1],digits=-1)
        ybreak2 = round(cibersortRange[2],digits=-1)
        
        p2 = ggplot(frame2, aes(y= CibersortPred ,x = RealCounts)) + geom_point(size = 6, color='red',shape=16) + geom_line(size=3,color='red',alpha=1/10) + theme_bw() %+replace% 
            theme(panel.background = element_rect(fill = NA)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                                      plot.title = element_text(size=35),
                                                                      axis.text = element_text(size=30)) + 
            ggtitle(paste0(names(estimates$estimates[i]),'\n(n genes = ',
                    estimates$rotations[[i]] %>% nrow, ')')) +
            scale_y_continuous(breaks = c(ybreak1, ybreak2), limits = c(cibersortRange[1], cibersortRange[2]+(cibersortRange[2] - cibersortRange[1])*.2))
        
        
        g1 <- ggplot_gtable(ggplot_build(p1))
        g2 <- ggplot_gtable(ggplot_build(p2))
        
        pp <- c(subset(g1$layout, name == "panel", se = t:r))
        g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                             pp$l, pp$b, pp$l)
        
        ia <- which(g2$layout$name == "axis-l")
        ga <- g2$grobs[[ia]]
        ax <- ga$children[[2]]
        ax$widths <- rev(ax$widths)
        ax$grobs <- rev(ax$grobs)
        ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
        g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
        g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
        png(paste0(outDir,names(estimates$estimates[i]),'.png'),width=600,height=600)
        grid.draw(g)
        dev.off()
    }
    correlations = t(as.data.frame(correlations))
    rownames(correlations) = names(estimates$estimates)
    #correlations[[1]] %<>% as.double
    #correlations[[2]] %<>% as.double
    correlations = data.frame(correlations, 
                              "n genes" = estimates$rotations %>% sapply(nrow),
                              check.names=FALSE)
    return(correlations)
}