# Sum of Inverse branch length

This statistic calculates per tip, the sum of the inverse of all branch
lengths of the shortest path between the root and the tip.

## Usage

``` r
inv_branch_dist(phy, add_one = FALSE)
```

## Arguments

- phy:

  phylo object or ltable

- add_one:

  should we take 1 / bl or 1 / (1 + bl) ? if TRUE, then 1 / (1 + bl) is
  used. Default: FALSE.

## References

Williams, P. H., Alonso-Alonso, P., Arbetman, M., Françoso, E.,
Ghisbain, G., Huang, J., Orr, M. C., Ren, Z.-X., Streinzer, M., T
hanoosing, C., Vandame, R., Waite, M., & Brace, S. (2026). Evolutionary
Tree for All Bumblebee Species World-Wide Estimated by Combining
Information from Fast-Evolving Genes, Slow-Evolving Genes, and Genomic
Data (Apidae, Bombus). Insects, 17(6), 540.
https://doi.org/10.3390/insects17060540
