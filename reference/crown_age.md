# Crown age of a tree.

In a reconstructed tree, obtaining the crown age is fairly
straightforward, and the function beautier::get_crown_age does a great
job at it. However, in a non-ultrametric tree, that function no longer
works. This function provides a functioning alternative.

## Usage

``` r
crown_age(phy)
```

## Arguments

- phy:

  phylo object or ltable

## Value

crown age
