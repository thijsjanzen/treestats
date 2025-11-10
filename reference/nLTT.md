# Normalized LTT statistic

The nLTT statistic calculates the sum of absolute differences in the
number of lineages over time, where both the number of lineages and the
time are normalized. The number of lineages is normalized by the number
of extant tips, whereas the time is normalized by the crown age. The
nLTT can only be calculated for reconstructed trees. Only use the
treestats version if you are very certain about the input data, and are
certain that performing nLTT is valid (e.g. your tree is ultrametric
etc). If you are less certain, use the nLTT function from the nLTT
package.

## Usage

``` r
nLTT(phy, ref_tree)
```

## Arguments

- phy:

  phylo object or ltable

- ref_tree:

  reference tree to compare with (should be same type as phy)

## Value

number of lineages

## References

Janzen, T., Höhna, S. and Etienne, R.S. (2015), Approximate Bayesian
Computation of diversification rates from molecular phylogenies:
introducing a new efficient summary statistic, the nLTT. Methods Ecol
Evol, 6: 566-575. https://doi.org/10.1111/2041-210X.12350

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
reference_tree <- ape::rphylo(n = 10, birth = 0.2, death = 0)
nLTT(simulated_tree, reference_tree)
#> [1] 0.09091369
nLTT(simulated_tree, simulated_tree) # should be zero.
#> [1] 0
```
