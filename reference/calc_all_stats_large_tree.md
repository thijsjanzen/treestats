# Apply all available tree statistics to a single tree, but exclude those statistics that require large amounts of memory.

this function applies all tree statistics available in this package to a
single tree, excluding statistics that require large amounts of memory,
because they are based on a distance matrix. The statistics EXCLUDED
are:

- laplacian spectrum

- variance of pairwise distance (vpd)

- maximum eigen vector value

- minimum eigenvalue of the Laplacian matrix

- maximum eigenvalue of the Laplacian matrix

- minimum eigenvalue of the adjacency matrix

- maximum eigenvalue of the adjacency matrix

## Usage

``` r
calc_all_stats_large_tree(phylo, normalize = FALSE)
```

## Arguments

- phylo:

  phylo object

- normalize:

  if set to TRUE, results are normalized (if possible) under either the
  Yule expectation (if available), or the number of tips

## Value

List with statistics
