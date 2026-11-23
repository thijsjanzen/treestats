# Stairs index

Calculates the staircase-ness measure, from the phyloTop package. The
staircase-ness reflects the number of subtrees that are imbalanced, e.g.
subtrees where the left child has more extant tips than the right child,
or vice versa.

## Usage

``` r
stairs(input_obj, normalization = "none")
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
