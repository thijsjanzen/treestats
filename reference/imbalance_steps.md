# Imbalance steps index

Calculates the number of moves required to transform the focal tree into
a fully imbalanced (caterpillar) tree. Higher value indicates a more
balanced tree.

## Usage

``` r
imbalance_steps(input_obj, normalization = FALSE)
```

## Arguments

- input_obj:

  phylo object or ltable

- normalization:

  if true, the number of steps taken is normalized by tree size, by
  dividing by the maximum number of moves required to move from a fully
  balanced to a fully imbalanced tree, which is \\N - log2(N) - 1\\,
  where N is the number of extant tips.

## Value

required number of moves
