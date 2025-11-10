# Height of a tree.

In a reconstructed tree, obtaining the tree height is fairly
straightforward, and the function beautier::get_crown_age does a great
job at it. However, in a non-ultrametric tree, that function no longer
works. Alternatively, taking the maximum value of
[adephylo::distRoot](https://rdrr.io/pkg/adephylo/man/distRoot.html)
will also yield the tree height (including the root branch), but will
typically perform many superfluous calculations and thus be slow.

## Usage

``` r
tree_height(phy)
```

## Arguments

- phy:

  phylo object

## Value

crown age
