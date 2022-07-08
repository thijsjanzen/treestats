context("max del width")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::max_del_width(focal_tree)
    a2 <- treebalance::maxDelW(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::max_del_width(focal_tree),
                           treestats::max_del_width(ltab))


    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    a1 <- treestats::max_del_width(focal_tree)
    a2 <- treebalance::maxDelW(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::max_del_width(focal_tree),
                           treestats::max_del_width(ltab))
  }
})
