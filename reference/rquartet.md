# Rquartet index.

The rquartet index counts the number of potential fully balanced rooted
subtrees of 4 tips in the tree. The function in treestats assumes a
bifurcating tree. For trees with polytomies, we use
treebalance::rquartedI, which can also take polytomies into account.

## Usage

``` r
rquartet(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  The index can be normalized by the expectation under the Yule ("yule")
  or PDA model ("pda").

## Value

rquartet index

## References

T. M. Coronado, A. Mir, F. Rosselló, and G. Valiente. A balance index
for phylogenetic trees based on rooted quartets. Journal of Mathematical
Biology, 79(3):1105-1148, 2019. doi: 10.1007/s00285-019-01377-w.
