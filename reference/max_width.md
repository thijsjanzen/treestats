# Maximum width of branch depths.

Calculates the maximum width, this is calculated by first collecting the
depth of each node and tip across the entire tree, where the depth
represents the distance (in nodes) to the root. Then, the width
represents the number of occurrences of each possible depth. The maximal
width then returns the maximum number of such occurences.

## Usage

``` r
max_width(phy, normalization = "none")
```

## Arguments

- phy:

  phylogeny or ltable

- normalization:

  "none" or "tips", in which case the resulting statistic is divided by
  the number of tips in the tree.

## Value

maximum width

## References

C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
transmission patterns. Evolution, Medicine, and Public Health,
2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
