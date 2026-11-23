# Getting started

## Using treestats

The treestats package provides an easy to use interface to calculate
summary statistics on phylogenetic trees. To obtain a list of all
supported summary statistics use:

``` r

list_statistics()
```

    ##  [1] "area_per_pair"          "average_leaf_depth"     "avg_ladder"            
    ##  [4] "avg_vert_depth"         "b1"                     "b2"                    
    ##  [7] "beta"                   "cherries"               "colless"               
    ## [10] "colless_corr"           "colless_quad"           "crown_age"             
    ## [13] "diameter"               "double_cherries"        "eigen_centrality"      
    ## [16] "eigen_centralityW"      "ew_colless"             "four_prong"            
    ## [19] "gamma"                  "i_stat"                 "il_number"             
    ## [22] "imbalance_steps"        "j_one"                  "j_stat"                
    ## [25] "laplace_spectrum_a"     "laplace_spectrum_e"     "laplace_spectrum_g"    
    ## [28] "laplace_spectrum_p"     "max_adj"                "max_betweenness"       
    ## [31] "max_clade_basis"        "max_closeness"          "max_closenessW"        
    ## [34] "max_del_width"          "max_depth"              "max_ladder"            
    ## [37] "max_laplace"            "max_width"              "mean_branch_length"    
    ## [40] "mean_branch_length_ext" "mean_branch_length_int" "mean_inv_branch_dist"  
    ## [43] "min_adj"                "min_laplace"            "mntd"                  
    ## [46] "mpd"                    "mw_over_md"             "nltt_base"             
    ## [49] "number_of_lineages"     "phylogenetic_div"       "pigot_rho"             
    ## [52] "pitchforks"             "psv"                    "rogers"                
    ## [55] "root_imbalance"         "rquartet"               "sackin"                
    ## [58] "sshape"                 "stairs"                 "stairs2"               
    ## [61] "symmetry_nodes"         "tot_coph"               "tot_internal_path"     
    ## [64] "tot_path"               "tree_height"            "treeness"              
    ## [67] "var_branch_length"      "var_branch_length_ext"  "var_branch_length_int" 
    ## [70] "var_depth"              "vpd"                    "wiener"

If your favourite summary statistic is missing, please let the
maintainer know, treestats is a dynamic package always under
development, and the maintainers are always looking for new statistics!

Given a phylogenetic tree, you can now use of the available functions to
calculate your summary statistic of choice. Let’s take for instance the
Colless statistic (and we generate a dummy tree):

``` r

phy <- ape::rphylo(n = 100, birth = 1, death = 0.1)

treestats::colless(phy)
```

    ## [1] 257

Looking at the documentation of the colless statistic
([`?colless`](https://thijsjanzen.github.io/treestats/reference/colless.md)),
we find that the function also includes options to normalize for size:
either ‘pda’ or ‘yule’:

``` r

treestats::colless(phy, normalization = "yule")
```

    ## [1] -0.9192387

### Multiple statistics

The treestats package supports calculating many statistics in one go.
For this, several functions have been set up aptly. Firstly, the
function `calc_all_stats` will calculate all statistics:

``` r

all_stats <- calc_all_stats(phy)
```

Similarly, we can also blanket apply all topology associated summary
statistics:

``` r

balance_stats <- calc_topology_stats(phy)
unlist(balance_stats)
```

    ##      area_per_pair average_leaf_depth         avg_ladder     avg_vert_depth 
    ##       1.255859e+01       7.630000e+00       2.250000e+00       6.673367e+00 
    ##                 b1                 b2               beta           cherries 
    ##       5.339845e+01       5.719727e+00       1.193623e+00       3.200000e+01 
    ##            colless       colless_corr       colless_quad           diameter 
    ##       2.570000e+02       5.297877e-02       4.235000e+03       2.000000e+01 
    ##    double_cherries   eigen_centrality         ew_colless         four_prong 
    ##       5.000000e+00       2.833317e-01       4.607934e-01       1.100000e+01 
    ##             i_stat          il_number    imbalance_steps              j_one 
    ##       4.712378e-01       3.600000e+01       8.700000e+01       8.707544e-01 
    ##    max_betweenness    max_clade_basis      max_closeness      max_del_width 
    ##       1.251500e+04       7.959908e+00                 NA       1.400000e+01 
    ##          max_depth         max_ladder          max_width         mw_over_md 
    ##       1.200000e+01       3.000000e+00       4.600000e+01       3.833333e+00 
    ##         pitchforks             rogers     root_imbalance           rquartet 
    ##       1.700000e+01       5.900000e+01       7.500000e-01       3.961488e+06 
    ##             sackin             sshape             stairs            stairs2 
    ##       7.630000e+02       1.589107e+02       5.959596e-01       6.650377e-01 
    ##     symmetry_nodes           tot_coph  tot_internal_path    tot_path_length 
    ##       5.900000e+01       6.686000e+03       5.650000e+02       1.328000e+03 
    ##          var_depth 
    ##       3.213100e+00
