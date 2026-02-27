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
  used. Default: TRUE

## Value

Distribution of the sum of inverse brnch length.
