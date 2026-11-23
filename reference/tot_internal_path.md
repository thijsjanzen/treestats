# Total internal path length

The total internal path length describes the sums of the depths of all
inner vertices of the tree.

## Usage

``` r
tot_internal_path(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "yule". When providing "yule", the statistic is normalized
  by the expectation under the Yule model, as given in Bauman &
  Rosenberg 2026.

## Value

Total internal path length

## References

Knuth, Donald E. The Art of Computer Programming: Fundamental
Algorithms, volume 1. Addison-Wesley Professional, 1997. Bauman, D. J.,
& Rosenberg, N. A. (2026). Basis statistics for phylogenetic tree shape:
tree balance formulas beyond expectations. arXiv preprint
arXiv:2610.40894.
