# \\J^1\\ index.

The \\J^1\\ index calculates the Shannon Entropy of a tree, where at
each node with two children, the Shannon Entropy is the sum of \\p_i
log_2(p_i)\\ over the two children \\i\\, and \\p_i\\ is \\L / (L +
R)\\, where L and R represent the number of tips connected to the two
daughter branches.

## Usage

``` r
j_one(input_obj)
```

## Arguments

- input_obj:

  phylo object or ltable

## Value

\\j^1\\ index

## References

Jeanne Lemant, Cécile Le Sueur, Veselin Manojlović, Robert Noble,
Robust, Universal Tree Balance Indices, Systematic Biology, Volume 71,
Issue 5, September 2022, Pages 1210–1224

https://doi.org/10.1093/sysbio/syac027
