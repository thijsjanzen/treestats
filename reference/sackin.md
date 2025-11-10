# Sackin index of (im)balance.

The Sackin index is calculated as the sum of ancestors for each of the
tips. Higher values indicate higher imbalance. Two normalizations are
available, where a correction is made for tree size, under either a Yule
expectation, or a pda expectation.

## Usage

``` r
sackin(phy, normalization = "none")
```

## Arguments

- phy:

  phylogeny or ltable

- normalization:

  normalization, either 'none' (default), "yule" or "pda".

## Value

Sackin index

## References

M. J. Sackin (1972). "Good" and "Bad" Phenograms. Systematic Biology.
21:225-226.

## Examples

``` r
simulated_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
sackin(balanced_tree)
#> [1] 34
sackin(unbalanced_tree) # should be much higher
#> [1] 54
```
