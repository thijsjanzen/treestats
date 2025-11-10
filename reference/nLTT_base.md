# Reference nLTT statistic

The base nLTT statistic can be used as a semi stand-alone statistic for
phylogenetic trees. However, please note that although this provides a
nice way of checking the power of the nLTT statistic without directly
comparing two trees, the nLTT_base statistic is not a substitute for
directly comparing two phylogenetic trees. E.g. one would perhaps
naively assume that \\nLTT(A, B) = \|nLTT(A, base) - nLTT(B, base)\\.
Indeed, in some cases this may hold true (when, for instance, all
normalized lineages of A are less than all normalized lineages of B),
but once the nLTT curve of A intersects the nLTT curve of B, this no
longer applies.

## Usage

``` r
nLTT_base(phy)
```

## Arguments

- phy:

  phylo object

## Value

number of lineages

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
nLTT_base(simulated_tree)
#> [1] 0.5417693
```
