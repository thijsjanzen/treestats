# ILnumber

The ILnumber is the number of internal nodes with a single tip child.
Higher values typically indicate a tree that is more unbalanced.

The ILnumber is the number of internal nodes with a single tip child, as
adapted from the phyloTop package.

## Usage

``` r
ILnumber(input_obj, normalization = "none")
```

## Arguments

- input_obj:

  phylo object or ltable

- normalization:

  "none" or "tips", in which case the result is normalized by dividing
  by N - 2, where N is the number of tips.

## Value

ILnumber
