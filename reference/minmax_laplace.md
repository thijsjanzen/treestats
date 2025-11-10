# Laplacian Matrix properties

Calculates the eigenvalues of the Laplacian Matrix, where the Laplacian
matrix is the matrix representation of a graph, in this case a
phylogeny. When the R package RSpectra is available, a faster
calculation can be used, which does not calculate all eigenvalues, but
only the maximum and minimum. As such, when using this option, the
vector of all eigenvalues is not returned.

## Usage

``` r
minmax_laplace(phy, use_rspectra = FALSE)
```

## Arguments

- phy:

  phylo object or ltable

- use_rspectra:

  boolean to indicate whether the helping package RSpectra should be
  used, in which case only the minimum and maximum values are returned

## Value

List with the minimum and maximum eigenvalues

## References

Chindelevitch, Leonid, et al. "Network science inspires novel tree shape
statistics." Plos one 16.12 (2021): e0259877.
