# Treeness statistic

Calculates the fraction of tree length on internal branches, also known
as treeness or stemminess

## Usage

``` r
treeness(phy)
```

## Arguments

- phy:

  phylo object or Ltable

## Value

sum of all internal branch lengths (e.g. branches not leading to a tip)
divided by the sum over all branch lengths.
