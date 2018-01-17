### Cross-laboratory analysis of brain cell type transcriptomes with applications to interpretation of bulk tissue data


B. Ogan Mancarci, Lilah Toker, Shreejoy Tripathy, Brenna Li, Brad Rocco, Etienne Sibille, Paul Pavlidis

This repository includes the code used for the NeuroExpresso paper. Installation and use requires devtools which can be acquired by doing

```r
install.packages('devtools')
```

The easiest way to repeat the analysis is to clone the library and running the analysis scripts from the root directory. From the root directory `devtools::load_all()` will make the functions and data available but is relatively slow. `devtools::install()` followed by `library(neuroExpressoAnalysis)` will allow faster access to the functions and data, this will also install dependencies.

