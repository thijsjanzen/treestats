# Total cophenetic index.

The total cophenetic index is the sum of the depth of the last common
ancestor of all pairs of leaves.

## Usage

``` r
tot_coph(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "yule", when "yule" is chosen, the statistic is divided by
  the Yule expectation

## Value

Total cophenetic index

## References

A. Mir, F. Rosselló, and L. Rotger. A new balance index for phylogenetic
trees. Mathematical Bio-sciences, 241(1):125-136, 2013. doi:
10.1016/j.mbs.2012.10.005.
