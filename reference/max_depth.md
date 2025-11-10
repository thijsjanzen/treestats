# Maximum depth metric

The maximum depth metric, measures the maximal path (in edges), between
the tips and the root.

## Usage

``` r
max_depth(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "tips", in which case the resulting statistic is divided by
  the number of tips in the tree.

## Value

Maximum depth (in number of edges)

## References

C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
transmission patterns. Evolution, Medicine, and Public Health,
2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
