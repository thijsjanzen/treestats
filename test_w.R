set.seed(5)
extant_tree <- ape::rphylo(n = 5, birth = 1, death = 0)

treestats::wiener(extant_tree, weight = TRUE)

calc_using_ape <- function(phy) {
  av2 <- ape::dist.nodes(phy)
  return(sum(av2[lower.tri(av2)]))
}

standard_test(treestats::wiener, calc_using_ape)

test_polytomies(treestats::wiener, calc_using_ape)

test_func <- function(phy) {
  treeCentrality::computeWienerIndex(phy, weight = TRUE)
}

a1 <- test_func(extant_tree)
a2 <- calc_using_ape(extant_tree)
a3 <- treestats::wiener(extant_tree, weight = TRUE)
cat(a1, a2, a3, "\n")

a1 <- test_func(complete_tree)
a2 <- calc_using_ape(complete_tree)
a3 <- treestats::wiener(complete_tree, weight = TRUE)
cat(a1, a2, a4, "\n")
# test_polytomies(test_func, calc_using_ape)





