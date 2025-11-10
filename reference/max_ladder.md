# Maximum ladder index

Calculate the maximum ladder index, from the phyloTop package. Higher
values indicate more unbalanced trees. To calculate the maximum ladder
index, first all potential ladders in the tree are calculated. A ladder
is defined as a sequence of nodes where one of the daughter branches is
a terminal branch, resulting in a 'ladder' like pattern. The maximum
ladder index then represents the longest ladder found among all observed
ladders in the tree.

## Usage

``` r
max_ladder(input_obj)
```

## Arguments

- input_obj:

  phylo object or ltable

## Value

longest ladder in the tree
