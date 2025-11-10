# s shape statistic of (im)balance.

The s shape statistic of imbalance (also known as the Blum statistic,
see [blum](https://thijsjanzen.github.io/treestats/reference/blum.md))
calculates the sum of \\log(N-1)\\ over all internal nodes, where N
represents the total number of extant tips connected to that node. An
alternative implementation can be found in the Castor R package.

## Usage

``` r
sshape(phy, normalization = FALSE)
```

## Arguments

- phy:

  phylogeny or ltable

- normalization:

  because the sshape statistic sums over all nodes, the resulting
  statistic tends to be correlated with the number of extant tips.
  Normalization can be performed by dividing by the number of extant
  tips.

## Value

s shape statistic of imbalance

## References

M. G. B. Blum and O. Francois (2006). Which random processes describe
the Tree of Life? A large-scale study of phylogenetic tree imbalance.
Systematic Biology. 55:685-691.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
  balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
  unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
  sshape(balanced_tree)
#> [1] 6.291569
  sshape(unbalanced_tree) # should be higher
#> [1] 12.80183
```
