# for gene permutations according to a certain tolerance level
# criteriaValue is a named vector
#' @export
selectRandom = function(gene,n , criteriaValue, tolerance = .05, invalids=c()){
    # invalids are there to prevent duplications
    criteriaValue %<>% .[!names(.) %in% invalids]
    medianExpPercentiles = ecdf(criteriaValue)(criteriaValue[gene])
    range = quantile(criteriaValue,
                     c(max(medianExpPercentiles-tolerance/2,0),    
                       min(medianExpPercentiles+tolerance/2,0.99)))
    
    eligible = which(criteriaValue>=range[1] & criteriaValue <= range[2])
    selection = names(criteriaValue[sample(eligible,n,replace=T)])
    return(selection)
}