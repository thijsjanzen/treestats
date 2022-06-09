context("tot_coph")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::tot_coph(focal_tree)
  a2 <- treebalance::totCophI(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::tot_coph(focal_tree),
                         treestats::tot_coph(ltab))




  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1 <- treestats::tot_coph(focal_tree)
  a2 <- treebalance::totCophI(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::tot_coph(focal_tree),
                         treestats::tot_coph(ltab))

})
