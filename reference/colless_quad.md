# Quadratic Colless index of (im)balance.

The Quadratic Colless index is calculated as the sum of \\(L - R)^2\\
over all nodes.

## Usage

``` r
colless_quad(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  A character string equals to "none" (default) for no normalization or
  "yule"

## Value

quadratic colless index

## References

Bartoszek, Krzysztof, et al. "Squaring within the Colless index yields a
better balance index." Mathematical Biosciences 331 (2021): 108503.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
colless_quad(balanced_tree)
#> [1] 6
colless_quad(unbalanced_tree) # should be higher
#> [1] 204
```
