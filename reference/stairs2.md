# Stairs2 index

Calculates the stairs2 measure, from the phyloTop package. The stairs2
reflects the imbalance at each node, where it represents the average
across measure at each node, the measure being \\min(L, R) / max(L,
R)\\, where L and R reflect the number of tips connected at the left (L)
and right (R) daughter.

## Usage

``` r
stairs2(input_obj, normalization = "none")
```

## Arguments

- input_obj:

  phylo object or ltable

- normalization:

  "none" or "yule". When providing "yule", the statistic is normalized
  by the expectation under the Yule model, as given in Bauman &
  Rosenberg 2026.

## Value

number of stairs

## References

Norström, Melissa M., et al. "Phylotempo: a set of r scripts for
assessing and visualizing temporal clustering in genealogies inferred
from serially sampled viral sequences." Evolutionary Bioinformatics 8
(2012): EBO-S9738. Bauman, D. J., & Rosenberg, N. A. (2026). Basis
statistics for phylogenetic tree shape: tree balance formulas beyond
expectations. arXiv preprint arXiv:2610.40894.
