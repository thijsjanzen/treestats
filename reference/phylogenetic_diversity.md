# Phylogenetic diversity at time point t

The phylogenetic diversity at time t is given by the total branch length
of the tree reconstructed up until time point t. Time is measured
increasingly, with the crown age equal to 0. Thus, the time at the
present is equal to the crown age.

## Usage

``` r
phylogenetic_diversity(input_obj, t = 0, extinct_tol = NULL)
```

## Arguments

- input_obj:

  phylo object or Ltable

- t:

  time point at which to measure phylogenetic diversity, alternatively a
  vector of time points can also be provided. Time is measured with 0
  being the present.

- extinct_tol:

  tolerance to determine if a lineage is extinct at time t. Default is
  1/100 \* smallest branch length of the tree.

## Value

phylogenetic diversity, or vector of phylogenetic diversity measures if
a vector of time points is used as input.

## References

Faith, Daniel P. "Conservation evaluation and phylogenetic diversity."
Biological conservation 61.1 (1992): 1-10.
