# Cherry index

Calculate the number of cherries, from the phyloTop package. A cherry is
a pair of sister tips.

## Usage

``` r
cherries(input_obj, normalization = "none")
```

## Arguments

- input_obj:

  phylo object or ltable

- normalization:

  "none", "yule", or "pda", the found number of cherries is divided by
  the expected number, following McKenzie & Steel 2000.

## Value

number of cherries

## References

McKenzie, Andy, and Mike Steel. "Distributions of cherries for two
models of trees." Mathematical biosciences 164.1 (2000): 81-92.
