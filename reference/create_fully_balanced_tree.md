# Create a fully balanced tree

This function takes an input phylogeny, and returns a phylogeny that is
most ideally balanced tree, whilst having the same branching times as
the original input tree. Please note that if the number of tips is not
even or not a power of two, the tree may not have perfect balance, but
the most ideal balance possible.

## Usage

``` r
create_fully_balanced_tree(phy)
```

## Arguments

- phy:

  phylo object

## Value

phylo phylo object

## Examples

``` r
phy <- ape::rphylo(n = 16, birth = 1, death = 0)
bal_tree <- treestats::create_fully_balanced_tree(phy)
treestats::colless(phy)
#> [1] 34
treestats::colless(bal_tree) # much lower
#> [1] 0
```
