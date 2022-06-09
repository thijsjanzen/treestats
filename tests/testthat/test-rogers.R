context("rogers")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  rog <- treestats::rogers(focal_tree)
  rog_check <- treebalance::rogersI(focal_tree)
  testthat::expect_equal(rog, rog_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::rogers(focal_tree),
                         treestats::rogers(ltab))


  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  rog <- treestats::rogers(focal_tree)
  rog_check <- treebalance::rogersI(focal_tree)
  testthat::expect_equal(rog, rog_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::rogers(focal_tree),
                         treestats::rogers(ltab))
})
