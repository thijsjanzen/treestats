# Variance of leaf depth statistic

The variance of leaf depth statistic returns the variance of depths
across all tips, where depth of a tip indicates the distance of the tip
to the root.

## Usage

``` r
var_leaf_depth(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "yule", when "yule" is chosen, the statistic is divided by
  the Yule expectation

## Value

Variance of leaf depths

## References

T. M. Coronado, A. Mir, F. Rosselló, and L. Rotger. On Sackin's original
proposal: the variance of the leaves' depths as a phylogenetic balance
index. BMC Bioinformatics, 21(1), 2020. doi: 10.1186/s12859-020-3405-1.
