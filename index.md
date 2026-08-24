# Treestats

![](https://github.com/thijsjanzen/treestats/blob/main/layout/hex_treestats.png?raw=true)

#### Measuring properties of phylogenetic trees

The **treestats** R package contains rapid, C++ based, functions to
calculate summary statistics on phylogenies. For some functions (but not
all, see below), the phylogenies are required to be ultrametric and/or
binary.

## Getting started

### Installation

To get started, you can either install from CRAN or use the latest
version from GitHub:

    install.packages("treestats") # install from CRAN

    # use the devtools package to install latest version from GitHub:
    install.packages("devtools")
    devtools::install_github("thijsjanzen/treestats")

### Basic usage

Given a tree (for example a simulated tree, as in the code example), you
can either access individual statistics, or calculate all currently
implemented statistics:

    focal_tree   <- ape::rphylo(n = 10, birth = 1, death = 0)
    colless_stat <- treestats::colless(focal_tree)
    all_stats    <- treestats::calc_all_stats(focal_tree)

## Rcpp

For all of the implemented statistics, the package provides Rcpp
versions that are much, much faster than their R sister functions.
Furthermore, some additional functions have been improved as well: \*
ape::branching.times \* DDD::phylo2L \* DDD::L2phylo

![](https://github.com/thijsjanzen/treestats/blob/cd3649833740eb7cdb23f722a2738cfd23bc4b10/layout/Figure_S3.png?raw=true)

### List of statistics

The following summary statistics are included:

|  |  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|----|
| Statistic | Information | Norma lization | Req Binary? | Req Ultrametric? | Req Rooted? | Sensitive to root position | Reference |
| area_per_pair | Topology | Yule | FALSE | FALSE | TRUE | YES | Lima et al., 2020 |
| average_leaf_depth | Topology | Yule | FALSE | FALSE | TRUE | YES | Shao & Sokal, 1990 |
| avg_ladder | Topology | None | TRUE | FALSE | TRUE | YES | Kendall et al., 2018 |
| avg_vert_depth | Topology | None | FALSE | FALSE | TRUE | YES | Colijn & Gardy, 2014 |
| b1 | Topology | Tips | FALSE | FALSE | TRUE | YES | Shao & Sokal, 1990 |
| b2 | Topology | Yule | FALSE | FALSE | TRUE | YES | Shao & Sokal, 1990 |
| beta | Topology | None | TRUE | FALSE | TRUE | YES | Aldous, 1996 |
| blum / sshape | Topology | None | FALSE | FALSE | TRUE | YES | Blum & François, 2006 |
| cherries | Topology | Yule | FALSE | FALSE | FALSE | YES | McKenzie et al., 1999 |
| colless | Topology | Yule | TRUE | FALSE | TRUE | YES | Colless, 1982 |
| colless_corr | Topology | None | TRUE | FALSE | TRUE | YES | Heard, 1992 |
| colless_quad | Topology | None | TRUE | FALSE | TRUE | YES | Bartoszek et al., 2021 |
| crown_age | Br times | None | FALSE | FALSE | TRUE | YES |  |
| diameter | Topology | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| double_cherries | Topology | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| eigen_centrality | Topology | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| eigen_centralityW | Topo + br len | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| ew_colless | Topology | None | TRUE | FALSE | TRUE | YES | Mooers & S. B. Heard, 1997 |
| four_prong | Topology | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| gamma | Br times | None | FALSE | TRUE | TRUE | YES | Pybus & Harvey, 2000 |
| i_stat | Topology | None | TRUE | FALSE | TRUE | YES | Fusco & Cronk, 1995 |
| il_number | Topology | Tips | FALSE | FALSE | TRUE | YES | Kendall et al., 2018 |
| imbalance_steps | Topology | Tips | TRUE | FALSE | TRUE | YES | Janzen & Etienne, 2024 |
| inv_branch_dist | Topo + br len | None | FALSE | TRUE | TRUE | YES | Williams et al. 2026 |
| j_one | Topology | None | TRUE | FALSE | TRUE | YES | Lemant et al., 2022 |
| j_stat | Topo + br len | None | FALSE | FALSE | FALSE | NO | Izsák & Papp, 2000 |
| laplace_spectrum_a | Topo + br len | None | FALSE | FALSE | FALSE | YES | Lewitus & Morlon, 2016 |
| laplace_spectrum_e | Topo + br len | None | FALSE | FALSE | FALSE | YES | Lewitus & Morlon, 2016 |
| laplace_spectrum_g | Topo + br len | None | FALSE | FALSE | FALSE | YES | Lewitus & Morlon, 2016 |
| laplace_spectrum_p | Topo + br len | None | FALSE | FALSE | FALSE | YES | Lewitus & Morlon, 2016 |
| max_adj | Topo + br len | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| max_betweenness | Topology | Tips | FALSE | FALSE | TRUE | YES | Chindelevitch et al., 2021 |
| max_closeness | Topology | Tips | FALSE | FALSE | FALSE | Small | Chindelevitch et al., 2021 |
| max_closenessW | Topo + br len | None | FALSE | FALSE | FALSE | Small | Chindelevitch et al., 2021 |
| max_del_width | Topology | Tips | FALSE | FALSE | TRUE | YES | Colijn & Gardy, 2014 |
| max_depth | Topology | Tips | FALSE | FALSE | TRUE | YES | Colijn & Gardy, 2014 |
| max_ladder | Topology | None | TRUE | FALSE | FALSE | YES | Kendall et al., 2018 |
| max_laplace | Topo + br len | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| max_width | Topology | Tips | FALSE | FALSE | TRUE | YES | Colijn & Gardy, 2014 |
| mean_branch_length | Topo + br len | None | FALSE | FALSE | FALSE | YES | Janzen & Etienne, 2017 |
| mean_branch_length_ext | Topo + br len | None | FALSE | FALSE | FALSE | NO | Saulnier et al., 2017 |
| mean_branch_length_int | Topo + br len | None | FALSE | FALSE | FALSE | YES | Saulnier et al., 2017 |
| min_adj | Topo + br len | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| min_laplace | Topo + br len | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |
| mntd | Topo + br len | None | FALSE | FALSE | FALSE | NO | Webb et al., 2002 |
| mpd | Topo + br len | Tips | FALSE | FALSE | FALSE | NO | Webb et al., 2002 |
| mw_over_md | Topology | None | FALSE | FALSE | TRUE | YES | Colijn & Gardy, 2014 |
| nltt_base | Br times | None | FALSE | TRUE | TRUE | YES | Janzen et al., 2015 |
| number_of_lineages | Topo + br len | None | FALSE | FALSE | FALSE | NO |  |
| phylogenetic_div | Topo + br len | None | FALSE | FALSE | TRUE | YES | Faith, 1992 |
| pigot_rho | Br times | None | FALSE | FALSE | TRUE | YES | Pigot et al., 2010 |
| pitchforks | Topology | Tips | FALSE | FALSE | FALSE | YES | Kendall et al., 2018 |
| psv | Topo + br len | Tips | FALSE | FALSE | FALSE | NO | Helmus et al., 2007 |
| rogers | Topology | Tips | TRUE | FALSE | FALSE | YES | Rogers, 1996 |
| root_imbalance | Topology | None | FALSE | FALSE | TRUE | YES | Guyer et al., 1993 |
| rquartet | Topology | Yule | FALSE | FALSE | TRUE | YES | Coronado et al., 2019 |
| sackin | Topology | Yule | FALSE | FALSE | TRUE | YES | Sackin, 1972 |
| stairs | Topology | None | TRUE | FALSE | TRUE | YES | Norström et al., 2012 |
| stairs2 | Topology | None | TRUE | FALSE | TRUE | YES | Norström et al., 2012 |
| symmetry_nodes | Topology | Tips | TRUE | FALSE | TRUE | YES | Kersting & Fischer, 2021 |
| tot_coph | Topology | Yule | FALSE | FALSE | TRUE | YES | Mir et al., 2013 |
| tot_internal_path | Topology | None | FALSE | FALSE | TRUE | YES | Knuth, 1997 |
| tot_path | Topology | None | FALSE | FALSE | TRUE | YES | Colijn & Gardy, 2014 |
| tree_height | Br times | None | FALSE | FALSE | TRUE | YES |  |
| treeness | Topo + br len | None | FALSE | FALSE | FALSE | YES | Astolfi & Zonta-Sgaramella, 1984 |
| var_branch_length | Topo + br len | None | FALSE | FALSE | FALSE | YES | Saulnier et al., 2017 |
| var_branch_length_ext | Topo + br len | None | FALSE | FALSE | FALSE | YES | Saulnier et al., 2017 |
| var_branch_length_int | Topo + br len | None | FALSE | FALSE | FALSE | YES | Saulnier et al., 2017 |
| var_depth | Topology | Yule | FALSE | FALSE | TRUE | YES | Coronado et al., 2020 |
| vpd | Topo + br len | None | FALSE | FALSE | FALSE | NO | Webb et al., 2002 |
| wiener | Topo + br len | None | FALSE | FALSE | FALSE | YES | Chindelevitch et al., 2021 |

## C++ Library

For the Rcpp improved summary statistics (excluding statistics that rely
on the calculation of eigen values, as these rely on the Rcpp
independent Eigen code), R independent C++ code is provided in the
inst/include folder. These can be independently linked by adding the
treestats package in the DESCRIPTION in both the LinkingTo and Depends
fields. Then, in your package, you can also calculate these functions.

Please note that for all functions, there are two versions available: 1)
based on input of a phylo object, which is typically one 2-column matrix
containing all edges, and a vector containing the edge lengths
(depending on which information is required to calculate the statistic).
2) based on input of an Ltable (Lineage table), which is a 4-column
matrix containing information on each species, being 1) birth time, 2)
parent species, 3) species label and 4) death time (or -1 if extant).

Ltable input can be useful when summary statistics are required for more
complicated simulation models.
