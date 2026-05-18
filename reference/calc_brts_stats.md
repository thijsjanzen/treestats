# Apply all tree statistics related to branching times to a single tree.

this function applies all tree statistics based on branching times to a
single tree (more or less ignoring topology), being:

- gamma

- pigot's rho

- mean branch length

- nLTT with empty tree

- treeness

- var branch length

- mean internal branch length

- mean external branch length

- var internal branch length

- var external branch length

## Usage

``` r
calc_brts_stats(phylo)
```

## Arguments

- phylo:

  phylo object

## Value

list with statistics
