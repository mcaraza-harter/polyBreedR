polyBreedR
================
Jeffrey Endelman

This R package was designed to facilitate the use of genome-wide markers
for breeding autotetraploid (4x) species, but it also functionality for
diploids. The software has been developed and tested using data from the
University of Wisconsin-Madison [potato breeding
program](http://potatobreeding.cals.wisc.edu).

To install and load the package:

``` r
install.packages("devtools")
library(devtools)
install_github("jendelman/polyBreedR", build_vignettes=FALSE)
library(polyBreedR)
```

### Vignettes

[Vignette
1](https://jendelman.github.io/polyBreedR/polyBreedR_Vignette1.html)
covers the making of tetraploid genotype calls and their utility for
pedigree curation. The pedigree functionality is based on [Endelman et
al. (2017)](https://doi.org/10.1007/s12230-016-9556-y). Financial
support has come from USDA NIFA Hatch Projects 1002731 and 1013047,
administered through the Wisconsin Agricultural Experiment Station at
UW-Madison.

[Vignette
2](https://jendelman.github.io/polyBreedR/polyBreedR_Vignette2.html)
covers the analysis of marker data for dihaploids, including checking
ploidy and parentage. Financial support has come from USDA NIFA Awards
2014-67013-22434 (AFRI) and 2019-51181-30021 (SCRI).

For a complete specification of the package functions, consult the
[reference
manual.](https://jendelman.github.io/polyBreedR/polyBreedR_Manual.pdf)
