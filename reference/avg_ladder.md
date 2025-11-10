# Average ladder index

Calculate the avgLadder index, from the phyloTop package. Higher values
indicate more unbalanced trees. To calculate the average ladder index,
first all potential ladders in the tree are calculated. A ladder is
defined as a sequence of nodes where one of the daughter branches is a
terminal branch, resulting in a 'ladder' like pattern. The average
ladder index then represents the average lenght across all observed
ladders in the tree.

## Usage

``` r
avg_ladder(input_obj)
```

## Arguments

- input_obj:

  phylo object or ltable

## Value

average number of ladders
