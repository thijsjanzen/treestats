# Maximum difference of widths of a phylogenetic tree

Calculates the maximum difference of widths of a phylogenetic tree.
First, the widths are calculated by collecting the depth of each node
and tip across the entire tree, where the depth represents the distance
(in nodes) to the root. Then, the width represents the number of
occurrences of each possible depth. Then, we take the difference between
each consecutive width, starting with the first width. The maximum
difference is then returned - whereas the original statistic designed by
Colijn and Gardy used the absolute maximum difference, we here use the
modified version as introduced in Fischer 2023: this returns the maximum
value, without absoluting negative widths. This ensures that this metric
is a proper (im)balance metric, following Fischer 2023.

## Usage

``` r
max_del_width(phy, normalization = "none")
```

## Arguments

- phy:

  phylogeny or ltable

- normalization:

  "none" or "tips", in which case the resulting statistic is divided by
  the number of tips in the tree.

## Value

maximum difference of widths

## References

C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease
transmission patterns. Evolution, Medicine, and Public Health,
2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.

Fischer, M., Herbst, L., Kersting, S., Kühn, A. L., & Wicke, K. (2023).
Tree Balance Indices: A Comprehensive Survey.
