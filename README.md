# treestats

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/treestats)](https://cran.r-project.org/package=treestats)
[![R-CMD-check](https://github.com/thijsjanzen/treestats/workflows/R-CMD-check/badge.svg)](https://github.com/thijsjanzen/treestats/actions)

Branch|CodeCov
---|---
master|[![codecov.io](https://codecov.io/gh/thijsjanzen/treestats/branch/master/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/treestats)
develop|[![codecov.io](https://codecov.io/gh/thijsjanzen/treestats/branch/develop/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/treestats)

## Description

R package that combines functions to calculate summary statistics on phylogenies.
The following summary statistics are included:
  - Aldous Beta
  - Sackin index
  - Colless index
  - Blum index
  - Crown age
  - Tree height
  - Pigot's Rho
  - Mean branch length
  - Number of lineages
  - Laplacian spectrum
  - nLTT
  - Gamma
  - Phylogenetic Diversity
  - avgLadder
  - ILnumber
  - pitchforks
  - stairs

## Rcpp
For some of these statistics, and some additional functions, the package provides Rcpp versions that 
are much, much faster than their R sister functions. These statistics are:
  - Aldous Beta
  - Sackin index
  - Colless index
  - nLTT
  - Gamma
  - Phylogenetic Diversity
  - RPANDA::spectR
  - ape::branching.times
  - DDD::phylo2L
  - DDD::L2phylo
  
## C++ Library
For the Rcpp improved summary statistics (excluding the RPANDA and DDD functions, 
as these are only partially captured in Rcpp), tree independent C++ code is provided 
in the inst/include folder. These can be independently linked by adding the treestats 
package in the DESCRIPTION in both the LinkingTo and Depends fields. Then, in your package,
you can also calculate these functions. 
Please note that for most functions, there are two versions available: 
1) based on input of a phylo object, which is typically one 2-column matrix containing all edges, and a vector containing the edge lengths.
2) based on input of an Ltable, which is a 4-column matrix containing information on each species, being 1) birth time, 2) parent species, 3) species label and 4) death time (or -1 if extant).

Ltable input can be useful when summary statistics are required for more complicated simulation models. 
