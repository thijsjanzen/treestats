# Stepwise increase the imbalance of a tree

the goal of this function is to increasingly imbalance a tree, by
changing the topology, one move at a time. It does so by re-attaching
terminal branches to the root lineage, through the ltable. In effect,
this causes the tree to become increasingly caterpillarlike. When
started with a balanced tree, this allows for exploring the gradient
between a fully balanced tree, and a fully unbalanced tree. Please note
that the algorithm will try to increase imbalance, until a fully
caterpillar like tree is reached, which may occur before unbal_steps is
reached. Three methods are available: "youngest", reattaches branches in
order of age, starting with the branch originating from the most recent
branching event and working itself through the tree. "Random" picks a
random branch to reattach. "Terminal" also picks a random branch, but
only from terminal branches (e.g. branches that don't have any daughter
lineages, which is maximized in a fully imbalanced tree).

## Usage

``` r
make_unbalanced_tree(
  init_tree,
  unbal_steps,
  group_method = "any",
  selection_method = "random"
)
```

## Arguments

- init_tree:

  starting tree to work with

- unbal_steps:

  number of imbalance generating steps

- group_method:

  choice of "any" and "terminal"

- selection_method:

  choice of "random", "youngest" and "oldest"

## Value

phylo object

## Examples

``` r
simulated_tree <- ape::rphylo(n = 16, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
intermediate_tree <- make_unbalanced_tree(balanced_tree, 8)
colless(balanced_tree)
#> [1] 0
colless(intermediate_tree) # should be intermediate value
#> [1] 71
colless(unbalanced_tree) # should be highest colless value
#> [1] 105
```
