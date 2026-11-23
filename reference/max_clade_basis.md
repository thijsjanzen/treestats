# Maximum Clade-basis statistic

The Maximum Clade-basis statistic measures consistency with the Yule
model, such that the statistic is maximized for a Yule tree. Values
lower than those expected under a Yule tree are indicative of the focal
tree not being the result of a Yule process. Although the statistic can
also be calculated for non-binary trees, the interpretation of the
statistic is unclear for non-binary trees.

## Usage

``` r
max_clade_basis(phy, normalization = "none")
```

## Arguments

- phy:

  phylogeny or ltable

- normalization:

  normalization, either 'none' (default) or yule

## Value

Statistic

## References

Bauman, D. J., & Rosenberg, N. A. (2026). Basis statistics for
phylogenetic tree shape: tree balance formulas beyond expectations.
arXiv preprint arXiv:2610.40894.
