# Corrected Colless index of (im)balance.

The Corrected Colless index is calculated as the sum of \\abs(L - R)\\
over all nodes, corrected for tree size by dividing over \\(n-1) \*
(n-2)\\, where n is the number of nodes.

## Usage

``` r
colless_corr(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  A character string equals to "none" (default) for no normalization or
  "yule", in which case the obtained index is divided by the Yule
  expectation.

## Value

corrected colless index

## References

Heard, Stephen B. "Patterns in tree balance among cladistic, phenetic,
and randomly generated phylogenetic trees." Evolution 46.6 (1992):
1818-1826.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
colless_corr(balanced_tree)
#> [1] 0.1111111
colless_corr(unbalanced_tree) # should be higher
#> [1] 1
```
