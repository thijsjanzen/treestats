# Wiener index

The Wiener index is defined as the sum of all shortest path lengths
between pairs of nodes in a tree.

## Usage

``` r
wiener(phy, normalization = FALSE, weight = TRUE)
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  if TRUE, the Wiener index is normalized by the number of nodes, e.g.
  by choose(n, 2), where n is the number of nodes.

- weight:

  if TRUE, branch lenghts are used.

## Value

Wiener index

## References

Chindelevitch, Leonid, et al. "Network science inspires novel tree shape
statistics." Plos one 16.12 (2021): e0259877. Mohar, B., Pisanski, T.
How to compute the Wiener index of a graph. J Math Chem 2, 267–277
(1988)
