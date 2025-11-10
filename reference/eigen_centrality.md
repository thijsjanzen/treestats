# Eigen vector centrality

Eigen vector centrality associates with each node \\v\\ the positive
value \\e(v)\\, such that: \\ \sum\_{e}^v w(uv) \* e(u) = \lambda \*
e(v) \\. Thus, \\e(v)\\ is the Perron-Frobenius eigenvector of the
adjacency matrix of the tree.

## Usage

``` r
eigen_centrality(phy, weight = TRUE, scale = FALSE, use_rspectra = FALSE)
```

## Arguments

- phy:

  phylo object or ltable

- weight:

  if TRUE, uses branch lengths.

- scale:

  if TRUE, the Eigenvector is rescaled

- use_rspectra:

  boolean to indicate whether the helping package RSpectra should be
  used, which is faster, but returns fewer eigen values.

## Value

List with the Eigen vector and the leading Eigen value

## References

Chindelevitch, Leonid, et al. "Network science inspires novel tree shape
statistics." Plos one 16.12 (2021): e0259877.
