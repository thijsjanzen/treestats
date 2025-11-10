# Gamma statistic

The gamma statistic measures the relative position of internal nodes
within a reconstructed phylogeny. Under the Yule process, the gamma
values of a reconstructed tree follow a standard normal distribution. If
gamma \> 0, the nodes are located more towards the tips of the tree, and
if gamma \< 0, the nodes are located more towards the root of the tree.
Only available for ultrametric trees.

## Usage

``` r
gamma_statistic(phy)
```

## Arguments

- phy:

  phylo object or ltable

## Value

gamma statistic

## References

Pybus, O. G. and Harvey, P. H. (2000) Testing macro-evolutionary models
using incomplete molecular phylogenies. Proceedings of the Royal Society
of London. Series B. Biological Sciences, 267, 2267–2272.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
gamma_statistic(simulated_tree) # should be around 0.
#> [1] -0.2218001
if (requireNamespace("DDD")) {
  ddd_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 7)$tes
  gamma_statistic(ddd_tree) # because of diversity dependence, should be < 0
}
#> Loading required namespace: DDD
#> [1] -3.672579
```
