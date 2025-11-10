# Rogers J index of (im)balance.

The Rogers index is calculated as the total number of internal nodes
that are unbalanced, e.g. for which both daughter nodes lead to a
different number of extant tips. in other words, the number of nodes
where \\L != R\\ (where L(R) is the number of extant tips of the Left
(Right) daughter node).

## Usage

``` r
rogers(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "tips", in which case the resulting statistic is divided by
  the number of tips - 2 (e.g. the maximum value of the rogers index for
  a tree).

## Value

Rogers index

## References

J. S. Rogers. Central Moments and Probability Distributions of Three
Measures of Phylogenetic Tree Imbalance. Systematic Biology,
45(1):99-110, 1996. doi: 10.1093/sysbio/45.1.99.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
rogers(balanced_tree)
#> [1] 3
rogers(unbalanced_tree) # should be higher
#> [1] 8
```
