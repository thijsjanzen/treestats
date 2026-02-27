# Adjancency Matrix properties

Calculates the eigenvalues of the Adjancency Matrix, where the Adjacency
matrix is a square matrix indicate whether pairs of vertices are
adjacent or not on a graph - here, entries in the matrix indicate
connections between nodes (and betweens nodes and tips). Entries in the
adjacency matrix are weighted by branch length. Then, using the
adjacency matrix, we calculate the spectral properties of the matrix,
e.g. the minimum positive and maximum positive eigenvalues of the
matrix. When the R package RSpectra is available, a faster calculation
can be used, which does not calculate all eigenvalues, but only the
maximum and minimum. As such, when using this option, the vector of all
eigenvalues is not returned.

## Usage

``` r
minmax_adj(phy, use_rspectra = TRUE, force_min = FALSE)
```

## Arguments

- phy:

  phylo object or ltable

- use_rspectra:

  boolean to indicate whether the helping package RSpectra should be
  used, in which case only the minimum and maximum values are returned

- force_min:

  boolean to indicate whether or not to also calculate the minimum
  positive eigenvalue if the tree is larger than 23170 tips. Default is
  FALSE, as this can be very slow, but in some cases might be of
  interest.

## Value

List with the minimum and maximum eigenvalues

## References

Chindelevitch, Leonid, et al. "Network science inspires novel tree shape
statistics." Plos one 16.12 (2021): e0259877.
