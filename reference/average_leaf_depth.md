# Average leaf depth statistic. The average leaf depth statistic is a normalized version of the Sackin index, normalized by the number of tips.

Average leaf depth statistic. The average leaf depth statistic is a
normalized version of the Sackin index, normalized by the number of
tips.

## Usage

``` r
average_leaf_depth(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "yule", in which case the statistic is divided by the
  expectation under the yule model, following Remark 1 in Coronado et
  al. 2020.

## Value

average leaf depth statistic

## References

M. Coronado, T., Mir, A., Rosselló, F. et al. On Sackin’s original
proposal: the variance of the leaves’ depths as a phylogenetic balance
index. BMC Bioinformatics 21, 154 (2020).
https://doi.org/10.1186/s12859-020-3405-1 K.-T. Shao and R. R. Sokal.
Tree balance. Systematic Zoology, 39(3):266, 1990. doi: 10.2307/2992186.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
average_leaf_depth(simulated_tree)
#> [1] 3.6
```
