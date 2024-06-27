# treestats

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/treestats)](https://cran.r-project.org/package=treestats)
[![](https://cranlogs.r-pkg.org/badges/grand-total/treestats)](https://cran.r-project.org/package=treestats)
[![](https://cranlogs.r-pkg.org/badges/treestats)](https://cran.r-project.org/package=treestats)
[![R-CMD-check](https://github.com/thijsjanzen/treestats/workflows/R-CMD-check/badge.svg)](https://github.com/thijsjanzen/treestats/actions)

Branch|CodeCov
---|---
master|[![codecov.io](https://codecov.io/gh/thijsjanzen/treestats/branch/master/graph/badge.svg)](https://app.codecov.io/gh/thijsjanzen/treestats)
develop|[![codecov.io](https://codecov.io/gh/thijsjanzen/treestats/branch/develop/graph/badge.svg)](https://app.codecov.io/gh/thijsjanzen/treestats)

## Description

R package that combines functions to calculate summary statistics on phylogenies.
The following summary statistics are included:
| Statistic              | Information               | Fischer   | Normalization | Assumes Ultrametric tree | Requires binary tree | Reference                        |
| ---------------------- | ------------------------- | --------- | ------------- | ------------------------ | -------------------- | -------------------------------- |
| area_per_pair          | Topology                  | Shape     | Yule          | NO                       | YES                  | Lima et al., 2020                |
| average_leaf_depth     | Topology                  | Imbalance | Yule          | NO                       | YES                  | Shao & Sokal, 1990               |
| avg_ladder             | Topology                  | Shape     | None          | NO                       | YES                  | Kendall et al., 2018             |
| avg_vert_depth         | Topology                  | Imbalance | None          | NO                       | NO                   | Colijn & Gardy, 2014             |
| b1                     | Topology                  | Balance   | Tips          | NO                       | NO                   | Shao & Sokal, 1990               |
| b2                     | Topology                  | Balance   | Yule          | NO                       | NO                   | Shao & Sokal, 1990               |
| beta                   | Topology                  | No index  | None          | NO                       | YES                  | Aldous, 1996                     |
| blum                   | Topology                  | Imbalance | None          | NO                       | YES                  | Blum & François, 2006            |
| cherries               | Topology                  | Shape     | Yule          | NO                       | YES                  | McKenzie et al., 1999            |
| colless                | Topology                  | Imbalance | Yule          | NO                       | YES                  | Colless, 1982                    |
| colless_corr           | Topology                  | Imbalance | None          | NO                       | YES                  | Heard, 1992                      |
| colless_quad           | Topology                  | Imbalance | None          | NO                       | YES                  | Bartoszek et al., 2021           |
| crown_age              | Branching times           | No index  | None          | NO                       | NO                   |                                  |
| diameter               | Topology                  | Shape     | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| double_cherries        | Topology                  | Shape     | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| eigen_centrality       | Topology                  | No index  | None          | NO                       | NO                   | Chindelevitch et al., 2021       |
| eigen_centralityW      | Topology + branch lengths | No index  | None          | NO                       | NO                   | Chindelevitch et al., 2021       |
| ew_colless             | Topology                  | Imbalance | None          | NO                       | YES                  | Mooers & S. B. Heard, 1997       |
| four_prong             | Topology                  | Shape     | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| gamma                  | Branching times           | No index  | None          | YES                      | NO                   | Pybus & Harvey, 2000             |
| i_stat                 | Topology                  | Imbalance | None          | NO                       | YES                  | Fusco & Cronk, 1995              |
| il_number              | Topology                  | Shape     | Tips          | NO                       | NO                   | Kendall et al., 2018             |
| imbalance_steps        | Topology                  | No index  | Tips          | NO                       | NO                   | Janzen & Etienne, 2024           |
| j_one                  | Topology                  | No index  | None          | NO                       | YES                  | Lemant et al., 2022              |
| j_stat                 | Topology + branch lengths | No index  | None          | NO                       | NO                   | Izsák & Papp, 2000               |
| laplace_spectrum_a     | Topology + branch lengths | No index  | None          | YES                      | NO                   | Lewitus & Morlon, 2016           |
| laplace_spectrum_e     | Topology + branch lengths | No index  | None          | YES                      | NO                   | Lewitus & Morlon, 2016           |
| laplace_spectrum_g     | Topology + branch lengths | No index  | None          | YES                      | NO                   | Lewitus & Morlon, 2016           |
| laplace_spectrum_p     | Topology + branch lengths | No index  | None          | YES                      | NO                   | Lewitus & Morlon, 2016           |
| max_adj                | Topology + branch lengths | No index  | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| max_betweenness        | Topology                  | Shape     | Tips          | NO                       | YES                  | Chindelevitch et al., 2021       |
| max_closeness          | Topology                  | Shape     | Tips          | NO                       | YES                  | Chindelevitch et al., 2021       |
| max_closenessW         | Topology + branch lengths | Shape     | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| max_del_width          | Topology                  | Shape     | Tips          | NO                       | NO                   | Colijn & Gardy, 2014             |
| max_depth              | Topology                  | Imbalance | Tips          | NO                       | NO                   | Colijn & Gardy, 2014             |
| max_ladder             | Topology                  | Shape     | None          | NO                       | YES                  | Kendall et al., 2018             |
| max_laplace            | Topology + branch lengths | No index  | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| max_width              | Topology                  | Balance   | Tips          | NO                       | NO                   | Colijn & Gardy, 2014             |
| mean_branch_length     | Topology + branch lengths | No index  | None          | NO                       | NO                   | Janzen & Etienne, 2017           |
| mean_branch_length_ext | Topology + branch lengths | No index  | None          | NO                       | NO                   | Saulnier et al., 2017            |
| mean_branch_length_int | Topology + branch lengths | No index  | None          | NO                       | NO                   | Saulnier et al., 2017            |
| min_adj                | Topology + branch lengths | No index  | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| min_laplace            | Topology + branch lengths | No index  | None          | NO                       | YES                  | Chindelevitch et al., 2021       |
| mntd                   | Topology + branch lengths | No index  | None          | NO                       | NO                   | Webb et al., 2002                |
| mpd                    | Topology + branch lengths | No index  | Tips          | NO                       | NO                   | Webb et al., 2002                |
| mw_over_md             | Topology                  | Balance   | None          | NO                       | NO                   | Colijn & Gardy, 2014             |
| nltt_base              | Branching times           | No index  | None          | YES                      | NO                   | Janzen et al., 2015              |
| number_of_lineages     | Topology + branch lengths | No index  | None          | NO                       | NO                   |                                  |
| phylogenetic_div       | Topology + branch lengths | No index  | None          | NO                       | NO                   | Faith, 1992                      |
| pigot_rho              | Branching times           | No index  | None          | YES                      | NO                   | Pigot et al., 2010               |
| pitchforks             | Topology                  | Shape     | Tips          | NO                       | NO                   | Kendall et al., 2018             |
| psv                    | Topology + branch lengths | No index  | Tips          | NO                       | NO                   | Helmus et al., 2007              |
| rogers                 | Topology                  | Imbalance | Tips          | NO                       | YES                  | Rogers, 1996                     |
| root_imbalance         | Topology                  | Shape     | None          | NO                       | YES                  | Guyer et al., 1993               |
| rquartet               | Topology                  | Balance   | Yule          | NO                       | NO                   | Coronado et al., 2019            |
| sackin                 | Topology                  | Imbalance | Yule          | NO                       | YES                  | Sackin, 1972                     |
| stairs                 | Topology                  | Imbalance | None          | NO                       | YES                  | Norström et al., 2012            |
| stairs2                | Topology                  | Balance   | None          | NO                       | YES                  | Norström et al., 2012            |
| symmetry_nodes         | Topology                  | Imbalance | Tips          | NO                       | YES                  | Kersting & Fischer, 2021         |
| tot_coph               | Topology                  | Imbalance | Yule          | NO                       | YES                  | Mir et al., 2013                 |
| tot_internal_path      | Topology                  | Imbalance | None          | NO                       | NO                   | Knuth, 1997                      |
| tot_path               | Topology                  | Imbalance | None          | NO                       | YES                  | Colijn & Gardy, 2014             |
| tree_height            | Branching times           | No index  | None          | NO                       | NO                   |                                  |
| treeness               | Topology + branch lengths | No index  | None          | NO                       | NO                   | Astolfi & Zonta-Sgaramella, 1984 |
| var_branch_length      | Topology + branch lengths | No index  | None          | NO                       | NO                   | Saulnier et al., 2017            |
| var_branch_length_ext  | Topology + branch lengths | No index  | None          | NO                       | NO                   | Saulnier et al., 2017            |
| var_branch_length_int  | Topology + branch lengths | No index  | None          | NO                       | NO                   | Saulnier et al., 2017            |
| var_depth              | Topology                  | Imbalance | Yule          | NO                       | NO                   | Coronado et al., 2020            |
| vpd                    | Topology + branch lengths | No index  | None          | NO                       | NO                   | Webb et al., 2002                |
| wiener                 | Topology + branch lengths | Shape     | None          | NO                       | YES                  | Chindelevitch et al., 2021       |

## Rcpp
For all of these statistics, the package provides Rcpp versions that 
are much, much faster than their R sister functions. Furthermore, some additional
functions have been improved as well:
  - ape::branching.times
  - DDD::phylo2L
  - DDD::L2phylo

[Figure_S3.pdf](https://github.com/user-attachments/files/16012805/Figure_S3.pdf)


## C++ Library
For the Rcpp improved summary statistics (excluding the RPANDA and DDD functions, 
as these are only partially captured in Rcpp), tree independent C++ code is provided 
in the inst/include folder. These can be independently linked by adding the treestats 
package in the DESCRIPTION in both the LinkingTo and Depends fields. Then, in your package,
you can also calculate these functions. 

Please note that for all functions, there are two versions available: 
1) based on input of a phylo object, which is typically one 2-column matrix containing all edges, and a vector containing the edge lengths (depending on which information is required to calculate the statistic).
2) based on input of an Ltable, which is a 4-column matrix containing information on each species, being 1) birth time, 2) parent species, 3) species label and 4) death time (or -1 if extant).

Ltable input can be useful when summary statistics are required for more complicated simulation models. 
