#' @export
translatePublishable = function(names){
    publishableNameDictionary$ShinyNames[
        match(names, 
              publishableNameDictionary$PyramidalDeep)]
}