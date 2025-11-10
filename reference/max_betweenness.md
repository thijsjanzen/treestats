# Maximum betweenness centrality.

Betweenness centrality associates with each node v, the two nodes u, w,
for which the shortest path between u and w runs through v, if the tree
were re-rooted at node v. Then, we report the node with maximum
betweenness centrality.

## Usage

``` r
max_betweenness(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "tips", if tips is chosen, the obtained maximum betweenness
  is normalized by the total amount of node pair combinations
  considered, e.g. (n-2)\*(n-1), where n is the number of tips.

## Value

Maximum Betweenness

## References

Chindelevitch, Leonid, et al. "Network science inspires novel tree shape
statistics." Plos one 16.12 (2021): e0259877.
