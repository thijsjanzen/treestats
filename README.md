[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/treestats)](https://cran.r-project.org/package=treestats)
[![](https://cranlogs.r-pkg.org/badges/grand-total/treestats)](https://cran.r-project.org/package=treestats)
[![](https://cranlogs.r-pkg.org/badges/treestats)](https://cran.r-project.org/package=treestats)
[![R-CMD-check](https://github.com/thijsjanzen/treestats/workflows/R-CMD-check/badge.svg)](https://github.com/thijsjanzen/treestats/actions)
[![codecov.io](https://codecov.io/gh/thijsjanzen/treestats/branch/main/graph/badge.svg)](https://app.codecov.io/gh/thijsjanzen/treestats) 


<img src="https://github.com/thijsjanzen/treestats/blob/main/layout/hex_treestats.png?raw=true" align="right" width="180" />

# Treestats 
#### Measuring properties of phylogenetic trees

The **treestats** R package contains rapid, C++ based, functions to
calculate summary statistics on phylogenies. For some functions (but not all, see below), the
phylogenies are required to be ultrametric and/or binary.

## Getting started

### Installation

To get started, you can either install from CRAN or use the latest
version from GitHub:

```         
install.packages("treestats") # install from CRAN

# use the devtools package to install latest version from GitHub:
install.packages("devtools")
devtools::install_github("thijsjanzen/treestats")
```

### Basic usage

Given a tree (for example a simulated tree, as in the code example), you
can either access individual statistics, or calculate all currently
implemented statistics:

```         
focal_tree   <- ape::rphylo(n = 10, birth = 1, death = 0)
colless_stat <- treestats::colless(focal_tree)
all_stats    <- treestats::calc_all_stats(focal_tree)
```

## Rcpp

For all of the implemented statistics, the package provides Rcpp versions that are
much, much faster than their R sister functions. Furthermore, some
additional functions have been improved as well: 
* ape::branching.times
* DDD::phylo2L
* DDD::L2phylo

![](https://github.com/thijsjanzen/treestats/blob/cd3649833740eb7cdb23f722a2738cfd23bc4b10/layout/Figure_S3.png?raw=true)

### List of statistics

The following summary statistics are included:
<table>
    <tr>
        <td>Statistic</td>
        <td>Information</td>
        <td>Normalization</td>
        <td>Fischer</td>
        <td>Req Binary?</td>
        <td>Req Ultrametric?</td>
        <td>Req Rooted?</td>
        <td>Sensitive to root position</td>
    </tr>
    <tr>
        <td>area_per_pair</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>average_leaf_depth</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>avg_ladder</td>
        <td>Topology</td>
        <td>None</td>
        <td>Shape</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>avg_vert_depth</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>b1</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Balance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>b2</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Balance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>beta</td>
        <td>Topology</td>
        <td>None</td>
        <td>No index</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>blum / sshape</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>cherries</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>colless</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>colless_corr</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>colless_quad</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>crown_age</td>
        <td>Br times</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>diameter</td>
        <td>Topology</td>
        <td>None</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>double_cherries</td>
        <td>Topology</td>
        <td>None</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>eigen_centrality</td>
        <td>Topology</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>eigen_centralityW</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>ew_colless</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>four_prong</td>
        <td>Topology</td>
        <td>None</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>gamma</td>
        <td>Br times</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>i_stat</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>il_number</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>imbalance_steps</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>No index</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>inv_branch_dist</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>j_one</td>
        <td>Topology</td>
        <td>None</td>
        <td>No index</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>j_stat</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>NO</td>
    </tr>
    <tr>
        <td>laplace_spectrum_a</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>laplace_spectrum_e</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>laplace_spectrum_g</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>laplace_spectrum_p</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>max_adj</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>max_betweenness</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>max_closeness</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>Small</td>
    </tr>
    <tr>
        <td>max_closenessW</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>Small</td>
    </tr>
    <tr>
        <td>max_del_width</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>max_depth</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>max_ladder</td>
        <td>Topology</td>
        <td>None</td>
        <td>Shape</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>max_laplace</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>max_width</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Balance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>mean_branch_length</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>mean_branch_length_ext</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>NO</td>
    </tr>
    <tr>
        <td>mean_branch_length_int</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>min_adj</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>min_laplace</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>mntd</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>NO</td>
    </tr>
    <tr>
        <td>mpd</td>
        <td>Topo + br len</td>
        <td>Tips</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>NO</td>
    </tr>
    <tr>
        <td>mw_over_md</td>
        <td>Topology</td>
        <td>None</td>
        <td>Balance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>nltt_base</td>
        <td>Br times</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>number_of_lineages</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>NO</td>
    </tr>
    <tr>
        <td>phylogenetic_div</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>pigot_rho</td>
        <td>Br times</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>pitchforks</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>psv</td>
        <td>Topo + br len</td>
        <td>Tips</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>NO</td>
    </tr>
    <tr>
        <td>rogers</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>root_imbalance</td>
        <td>Topology</td>
        <td>None</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>rquartet</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Balance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>sackin</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>stairs</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>stairs2</td>
        <td>Topology</td>
        <td>None</td>
        <td>Balance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>symmetry_nodes</td>
        <td>Topology</td>
        <td>Tips</td>
        <td>Imbalance</td>
        <td>TRUE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>tot_coph</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>tot_internal_path</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>tot_path</td>
        <td>Topology</td>
        <td>None</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>tree_height</td>
        <td>Br times</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>treeness</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>var_branch_length</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>var_branch_length_ext</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>var_branch_length_int</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>var_depth</td>
        <td>Topology</td>
        <td>Yule</td>
        <td>Imbalance</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>TRUE</td>
        <td>YES</td>
    </tr>
    <tr>
        <td>vpd</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>No index</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>NO</td>
    </tr>
    <tr>
        <td>wiener</td>
        <td>Topo + br len</td>
        <td>None</td>
        <td>Shape</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>FALSE</td>
        <td>YES</td>
    </tr>
</table>

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
