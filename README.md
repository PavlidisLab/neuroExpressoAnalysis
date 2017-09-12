### NeuroExpresso: A cross-laboratory database of brain cell-type expression profiles with applications to marker gene identification and bulk brain tissue transcriptome interpretation


B. Ogan Mancarci, Lilah Toker, Shreejoy Tripathy, Brenna Li, Brad Rocco, Etienne Sibille, Paul Pavlidis

This repository includes the code used for the NeuroExpresso paper. Installation and use requires devtools which can be acquired by doing

```r
install.packages('devtools')

```

The package itself depends on several personal libraries. Install them by
```r
devtools::install_github('oganm/ogbox')
devtools::install_github('oganm/homologene')
devtools::install_github('oganm/geneSynonym')
devtools::install_github('oganm/allenBrain')
devtools::install_github('oganm/markerGeneProfile')
```

There are also bioconductor packages in the dependency list that can be acquired by
```r
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("preprocessCore")
biocLite("limma")
biocLite("oligo")
biocLite("sva")
```

The easiest way to repeat the analysis is to clone the library and running the analysis scripts from the root directory. From the root directory `devtools::load_all()` will make the functions and data available but is relatively slow. `devtools::install()` followed by `library(neuroExpressoAnalysis)` will allow faster access to the functions and data.