context("mean I")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::mean_i(focal_tree)
  a2 <- treebalance::IbasedI(focal_tree, method = "mean",
                             correction = "prime", logs = FALSE)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::mean_i(focal_tree),
                         treestats::mean_i(ltab))


  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1 <- treestats::mean_i(focal_tree)
  a2 <- treebalance::IbasedI(focal_tree, method = "mean",
                             correction = "prime", logs = FALSE)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::mean_i(focal_tree),
                         treestats::mean_i(ltab))
})
