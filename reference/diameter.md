# Diameter statistic

The Diameter of a tree is defined as the maximum length of a shortest
path. When taking branch lengths into account, this is equal to twice
the crown age. When the tree is unrooted, we add 1 to the unweighted
diameter, to reflect traversing the (virtual) root.

## Usage

``` r
diameter(phy, weight = FALSE)
```

## Arguments

- phy:

  phylo object or ltable

- weight:

  if TRUE, uses branch lengths.

## Value

Diameter

## References

Chindelevitch, Leonid, et al. "Network science inspires novel tree shape
statistics." PloS one 16.12 (2021): e0259877.
