# Mean Pairwise distance

Fast function using C++ to calculate the mean pairwise distance, using
the fast algorithm by Constantinos, Sandel & Cheliotis (2012).

## Usage

``` r
mean_pair_dist(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "tips", in which case the obtained mean pairwise distance is
  normalized by the factor 2log(n), where n is the number of tips.

## Value

Mean pairwise distance

## References

Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002. Phylogenies and
community ecology. Annual Review of Ecology and Systematics 33:475-505.

Tsirogiannis, Constantinos, Brody Sandel, and Dimitris Cheliotis.
"Efficient computation of popular phylogenetic tree measures."
Algorithms in Bioinformatics: 12th International Workshop, WABI 2012,
Ljubljana, Slovenia, September 10-12, 2012. Proceedings 12. Springer
Berlin Heidelberg, 2012.
