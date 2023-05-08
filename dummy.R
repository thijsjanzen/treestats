library(treestats)

set.seed(42)

## with extinct lineages:
focal_tree <- ape::rphylo(n = 100, birth = 0.9, death = 0.2, fossils = TRUE)

colless <- treestats::colless(focal_tree)
colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree))

print(colless)
print(colless_check)
