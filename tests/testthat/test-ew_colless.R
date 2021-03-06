context("ew_colless")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    colless <- treestats::ew_colless(focal_tree)
    colless_check <- treebalance::ewCollessI(focal_tree)
    testthat::expect_equal(colless, colless_check)

    focal_ltab <- treestats::phylo_to_l(focal_tree)

    colless <- treestats::ew_colless(focal_ltab)
    testthat::expect_equal(colless, colless_check)

    ## with extinct lineages:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    colless <- treestats::ew_colless(focal_tree)
    colless_check <- treebalance::ewCollessI(focal_tree)
    testthat::expect_equal(colless, colless_check)

    focal_ltab <- treestats::phylo_to_l(focal_tree)

    colless <- treestats::ew_colless(focal_ltab)
    testthat::expect_equal(colless, colless_check)
  }
})
