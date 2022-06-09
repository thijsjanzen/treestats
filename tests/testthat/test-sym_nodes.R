context("symNodes")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::sym_nodes(focal_tree)
  a1_check <- treebalance::symNodesI(focal_tree)
  testthat::expect_equal(a1, a1_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::sym_nodes(focal_tree),
                         treestats::sym_nodes(ltab))


  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1 <- treestats::sym_nodes(focal_tree)
  a1_check <- treebalance::symNodesI(focal_tree)
  testthat::expect_equal(a1, a1_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::sym_nodes(focal_tree),
                         treestats::sym_nodes(ltab))
})
