# Create an unbalanced tree (caterpillar tree)

This function takes an input phylogeny, and returns a phylogeny that is
a perfectly imbalanced tree (e.g. a full caterpillar tree), that has the
same branching times as the original input tree.

## Usage

``` r
create_fully_unbalanced_tree(phy)
```

## Arguments

- phy:

  phylo object

## Value

phylo phylo object

## Examples

``` r
phy <- ape::rphylo(n = 16, birth = 1, death = 0)
bal_tree <- treestats::create_fully_unbalanced_tree(phy)
treestats::colless(phy)
#> [1] 31
treestats::colless(bal_tree) # much higher
#> [1] 105
```
