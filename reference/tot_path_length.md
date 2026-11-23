# Total path length

The total path length describes the sums of the depths of all vertices
of the tree. This can also be calculated as the sum of the total
internal path and the sackin index.

## Usage

``` r
tot_path_length(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "yule". When providing "yule", the statistic is normalized
  by the expectation under the Yule model, as given in Bauman &
  Rosenberg 2026.

## Value

Total path length

## References

C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
transmission patterns. Evolution, Medicine, and Public Health,
2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018. Bauman,
D. J., & Rosenberg, N. A. (2026). Basis statistics for phylogenetic tree
shape: tree balance formulas beyond expectations. arXiv preprint
