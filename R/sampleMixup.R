sampleMixup = function(expression, groups, p.value =  1e-5, lfc = log(5,base=2)){
    
    assertthat::assert_that(ncol(expression) == length(groups))
    uniGroups = unique(groups)
    
    pairwise = combn(uniGroups,2)
    
    difs = vector(mode = 'list', length = ncol(pairwise))
    
    for (i in 1:ncol(pairwise)){
        subsetExpr = expr[,groups %in% pairwise[,i]]
        subsetExpr = as.matrix(subsetExpr)
        # rownames(subsetExpr) = genes$Gene.Symbol
        
        subsetGroups = groups[groups %in% pairwise[,i]]

        mm = model.matrix(~ subsetGroups,data.frame(subsetGroups))
        fit <- lmFit(subsetExpr, mm)
        fit <- eBayes(fit)
        dif = topTable(fit, coef=colnames(fit$design)[2],
                       lfc = lfc,
                       number = Inf, 
                       p.value = p.value)
        if (nrow(dif)==0){
            print(paste('no genes found for',paste(pairwise[,i],collapse='/') ,'pair'))
            next
        }
        difs[[i]] = rownames(dif)
        #difs[[i]] = dif$ID
        print(i)
    }
    
    outliers = lapply(1:len(uniGroups), function(i){
        print(i)
        relevant = unlist(difs[apply(pairwise,2,function(x){
            uniGroups[i] %in%  x  
        })])
        subsetExpr = expression[,groups %in% uniGroups[i]]
        subsetExpr = subsetExpr[rownames(subsetExpr) %in% relevant,]
        pca = prcomp(t(subsetExpr))
        names(boxplot.stats(pca$x[,1])$out)
    })
    
    return(unlist(outliers))
}