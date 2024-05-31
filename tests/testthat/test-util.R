context("util-functions")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 5, birth = 1, death = 0)

  c1 <- treestats::colless(focal_tree, "yule")
  c2 <- treestats::colless(focal_tree, "Yule")
  testthat::expect_equal(c1, c2)

  c1 <- treestats::colless(focal_tree, "pda")
  c2 <- treestats::colless(focal_tree, "PDA")
  testthat::expect_equal(c1, c2)

  sn1 <- treestats::sym_nodes(focal_tree, "tips")
  sn2 <- treestats::sym_nodes(focal_tree, "Tips")
  testthat::expect_equal(sn1, sn2)
})
