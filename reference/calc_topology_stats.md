# Calculate all topology based statistics for a single tree

this function calculates all tree statistics based on topology available
in this package for a single tree, being:

- area_per_pair

- average_leaf_depth

- avg_ladder

- avg_vert_depth

- b1

- b2

- beta

- blum

- cherries

- colless

- colless_corr

- colless_quad

- diameter

- double_cherries

- eigen_centrality

- ew_colless

- four_prong

- i_stat

- il_number

- imbalance_steps

- j_one

- max_betweenness

- max_closeness

- max_del_width

- max_depth

- max_ladder

- max_width

- mw_over_md

- pitchforks

- rogers

- root_imbalance

- rquartet

- sackin

- stairs

- stairs2

- symmetry_nodes

- tot_coph

- tot_internal_path

- tot_path_length

- var_depth

## Usage

``` r
calc_topology_stats(phylo, normalize = FALSE)
```

## Arguments

- phylo:

  phylo object

- normalize:

  if set to TRUE, results are normalized (if possible) under either the
  Yule expectation (if available), or the number of tips

## Value

list with statistics
