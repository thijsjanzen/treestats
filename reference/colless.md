# Colless index of (im)balance.

The Colless index is calculated as the sum of \\abs(L - R)\\ over all
nodes, where L (or R) is the number of extant tips associated with the L
(or R) daughter branch at that node. Higher values indicate higher
imbalance. Two normalizations are available, where a correction is made
for tree size, under either a yule expectation, or a pda expectation.

## Usage

``` r
colless(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  A character string equals to "none" (default) for no normalization or
  one of "pda" or "yule".

## Value

colless index

## References

Colless D H. 1982. Review of: Phylogenetics: The Theory and Practice of
Phylogenetic Systematics. Systematic Zoology 31:100-104.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
colless(balanced_tree)
#> [1] 4
colless(unbalanced_tree) # should be higher
#> [1] 36
```
