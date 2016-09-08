library(viridis)

#' Return default colors for the cell types
#' @return A named vector. Due to existance of different naming schemes, 
#' same cell type can be referred to multiple times under different names.
#' @export
cellColors = function(){
    
    coloring = c(Oligo = 'darkgreen',
                 Oligodendrocyte = 'darkgreen',
                 Bergmann = 'palegreen',
                 MotorCholin = 'darkorange4',
                 SpinalChordCholinergic = 'darkorange4',
                 BrainstemCholin = 'darkorange4',
                 'Midbrain Cholinergic' = 'darkorange4',
                 ForebrainCholin = 'darkorange',
                 'Forebrain Cholinergic' = 'darkorange',
                 ThalamusCholin = 'darkorange4',
                 'Thalamus Cholinergic' = 'darkorange4',
                 Cholin = 'darkorange',
                 Spiny = 'blanchedalmond',
                 Gluta = 'slategray',
                 Basket = 'mediumpurple4',
                 Golgi = 'orchid',
                 Pyramidal = 'turquoise',
                 Purkinje = 'purple',
                 Inter = 'pink',
                 CerebGranule = 'thistle',
                 DentateGranule = 'thistle3',
                 Microglia = 'cornsilk4',
                 Microglia_activation = 'cornsilk3',
                 Microglia_deactivation = 'cornsilk2',
                 # Gaba = 'firebrick4',
                 Astrocyte = 'darkgoldenrod1',
                 GabaPV = 'firebrick2',
                 'FS Basket (G42)' = 'firebrick2',
                 Stem = 'blue' ,
                 Ependymal = 'orange',
                 Serotonergic = 'darkolivegreen',
                 Hypocretinergic = 'cadetblue',
                 Dopaminergic = 'gray0',
                 Th_positive_LC = 'blueviolet',
                 Noradrenergic = 'blueviolet',
                 GabaVIPReln = 'firebrick4',
                 'VIPReln (G30)' = 'firebrick4',
                 GabaRelnCalb = 'firebrick3',
                 'Martinotti (GIN)' = 'firebrick3',
                 GabaSSTReln = 'firebrick1',
                 GabaReln = 'firebrick',
                 GabaOxtr = 'firebrick2',
                 GabaHtr3a = 'darkred',
                 GabaVIPReln = 'firebrick4',
                 GabaReln = 'firebrick',
                 Pyramidal_Thy1 = 'turquoise',
                 PyramidalCorticoThalam = 'blue',
                 Pyramidal_Glt_25d2 = 'blue4',
                 Pyramidal_S100a10 ='deepskyblue3',
                 Glia = viridis(2)[1],
                 Neuron = viridis(2)[2]
                 
    )
    
    return(coloring)
}