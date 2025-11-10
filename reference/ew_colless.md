# Equal weights Colless index of (im)balance.

The equal weights Colless index is calculated as the sum of \\abs(L - R)
/ (L + R - 2)\\ over all nodes where \\L + R \> 2\\, where L (or R) is
the number of extant tips associated with the L (or R) daughter branch
at that node. Maximal imbalance is associated with a value of 1.0. The
ew_colless index is not sensitive to tree size.

## Usage

``` r
ew_colless(phy)
```

## Arguments

- phy:

  phylo object or ltable

## Value

colless index

## References

A. O. Mooers and S. B. Heard. Inferring Evolutionary Process from
Phylogenetic Tree Shape. The Quarterly Review of Biology, 72(1), 1997.
doi: 10.1086/419657.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
ew_colless(balanced_tree)
#> [1] 0.28125
ew_colless(unbalanced_tree) # should be higher
#> [1] 1
```
