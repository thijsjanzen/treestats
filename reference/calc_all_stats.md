# Apply all available tree statistics to a single tree

this function applies all tree statistics available in this package to a
single tree, being:

- gamma

- Sackin

- Colless

- corrected Colless

- quadratic Colless

- Aldous' beta statistic

- Blum

- crown age

- tree height

- Pigot's rho

- number of lineages

- nLTT with empty tree

- phylogenetic diversity

- avgLadder index

- cherries

- double cherries

- ILnumber

- pitchforks

- stairs

- stairs2

- laplacian spectrum

- B1

- B2

- area per pair (aPP)

- average leaf depth (aLD)

- I statistic

- ewColless

- max Delta Width (maxDelW)

- maximum of Depth

- variance of Depth

- maximum Width

- Rogers

- total Cophenetic distance

- symmetry Nodes

- mean of pairwise distance (mpd)

- variance of pairwise distance (vpd)

- Phylogenetic Species Variability (psv)

- mean nearest taxon distance (mntd)

- J statistic of entropy

- rquartet index

- Wiener index

- max betweenness

- max closeness

- diameter, without branch lenghts

- maximum eigen vector value

- mean branch length

- variance of branch length

- mean external branch length

- variance of external branch length

- mean internal branch length

- variance of internal branch length

- number of imbalancing steps

- j_one statistic

- treeness statistic

For the Laplacian spectrum properties, four properties of the eigenvalue
distribution are returned: 1) asymmetry, 2) peakedness, 3) log(principal
eigenvalue) and 4) eigengap. Please notice that for some very small or
very large trees, some of the statistics can not be calculated. The
function will report an NA for this statistic, but will not break, to
facilitate batch analysis of large numbers of trees.

## Usage

``` r
calc_all_stats(phylo, normalize = FALSE)
```

## Arguments

- phylo:

  phylo object

- normalize:

  if set to TRUE, results are normalized (if possible) under either the
  Yule expectation (if available), or the number of tips

## Value

List with statistics
