# Blum index of (im)balance.

The Blum index of imbalance (also known as the s-shape statistic (see
[sshape](https://thijsjanzen.github.io/treestats/reference/sshape.md)))
calculates the sum of \\log(N-1)\\ over all internal nodes, where N
represents the total number of extant tips connected to that node. An
alternative implementation can be found in the Castor R package.

## Usage

``` r
blum(phy, normalization = "none", log_base = 2)
```

## Arguments

- phy:

  phylogeny or ltable

- normalization:

  "none" or "yule", in which case the acquired result is divided by the
  expectation for the Yule model.

- log_base:

  the logarithm base that will be used. For binary trees it is common to
  use base 2.

## Value

Blum index of imbalance

## References

M. G. B. Blum and O. Francois (2006). Which random processes describe
the Tree of Life? A large-scale study of phylogenetic tree imbalance.
Systematic Biology. 55:685-691.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
  balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
  unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
  blum(balanced_tree)
#> [1] 9.076816
  blum(unbalanced_tree) # should be higher
#> [1] 18.46913
```
