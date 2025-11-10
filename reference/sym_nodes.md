# Symmetry nodes metric

Balance metric that returns the total number of internal nodes that are
not-symmetric (confusingly enough). A node is considered symmetric when
both daughter trees have the same topology, measured as having the same
sum of depths, where depth is measured as the distance from the root to
the node/tip.

## Usage

``` r
sym_nodes(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "tips", in which case the resulting statistic is divided by
  the number of tips - 2 (e.g. the maximum value of the symmetry nodes
  index for a tree).

## Value

Maximum depth (in number of edges)

## References

S. J. Kersting and M. Fischer. Measuring tree balance using symmetry
nodes — A new balance index and its extremal properties. Mathematical
Biosciences, page 108690, 2021. ISSN 0025-5564.
doi:https://doi.org/10.1016/j.mbs.2021.108690
