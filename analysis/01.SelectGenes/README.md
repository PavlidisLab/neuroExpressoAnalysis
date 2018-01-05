Marker gene changelog
======================

### Marker Genes

* [Markers_1.1](Markers_1.1/README.md)
* [Markers_1.0](Markers_1.0/README.md)

### Markers_1.1

This update introduces minor fixes to the gene list due to the following factors

* A bug in the function that places samples into brain region. This resulted in oligodendrocytes
isolated from the mouse neocortex to be added while selecting genes for cerebellar cell types even
though we intended them to be not included since we had cerebellar oligodendrocyte samples.


* An update to the sample metada that changed which probes are selected to 
represent a gene during normalization. Probes selected for a 116 genes was effected
due to this change. 16 new genes were included in the analysis due to the newly selected probe
being above out selection threshold.

For a complete list of genes added or removed from the "Combined" gene list see [this file](Markers_1.1/comparison_CellTypes-Markers_1.0-Markers_1.1.md).

Below is the number of genes added or removed from all cell types across brain regions

|                       | added| removed| Total genes in final list|
|:----------------------|-----:|-------:|-------------------------:|
|Astrocyte              |    13|      19|                       426|
|Microglia              |    37|      25|                       844|
|Microglia_activation   |     6|       5|                       187|
|Microglia_deactivation |    12|       8|                       279|
|Noradrenergic          |     2|       4|                       135|
|Oligo                  |    18|       9|                       435|
|ForebrainCholin        |     4|       0|                       106|
|BrainstemCholin        |     1|       1|                        40|
|Basket                 |     1|       1|                         7|
|Bergmann               |     8|       2|                        71|
|CerebGranule           |     0|       5|                        12|
|Golgi                  |     4|       1|                        30|
|Purkinje               |     2|       1|                        48|
|DentateGranule         |     1|       0|                        22|
|GabaPV                 |     0|       1|                        30|
|GabaRelnCalb           |     1|       0|                        24|
|GabaVIPReln            |     1|       0|                        49|
|GabaSSTReln            |     1|       0|                        55|
|SpinalCordCholinergic  |     1|       0|                       128|
|Spiny                  |     1|       2|                        74|
|Ependymal              |     1|       0|                        68|
|Dopaminergic           |     1|       2|                        60|




### Markers_1.0

Initial version submitted with paper