#' Quick fix for ugly data tables
#' @export
fixTable = function(x){
    x %>%  apply(2,function(x){x %>% sapply(function(y){
        if(is.na(y)){
            return('genes<2')
        } else if(is.numeric(y)){
            y %<>% round(digits=3)
            if(y==0){
                return('p<0.001')
            }
            return(y %>% round(digits=3))
        }
        return(y)
    })})
}