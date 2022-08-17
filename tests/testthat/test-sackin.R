context("sackin")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  sackin <- treestats::sackin(focal_tree)
  sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(focal_tree))
  testthat::expect_equal(sackin, sackin_check)

  sackin <- treestats::sackin(focal_tree, normalization = "yule")
  sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(focal_tree),
                                      norm = "yule")
  testthat::expect_equal(sackin, sackin_check, tol = 1e-5)

  sackin <- treestats::sackin(focal_tree, normalization = "pda")
  sackin_check <- apTreeshape::sackin(apTreeshape::as.treeshape(focal_tree),
                                      norm = "pda")
  testthat::expect_equal(sackin, sackin_check, tol = 1e-5)


  # test ltable functionality
  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::sackin(focal_tree),
                        treestats::sackin(ltab))

  testthat::expect_equal(treestats::sackin(focal_tree, normalization = "yule"),
                         treestats::sackin(ltab, normalization = "yule"))

  testthat::expect_equal(treestats::sackin(focal_tree, normalization = "pda"),
                         treestats::sackin(ltab, normalization = "pda"))

  if (requireNamespace("treebalance")) {
    sackin <- treestats::sackin(focal_tree)
    sackin_check <- treebalance::sackinI(focal_tree)
    testthat::expect_equal(sackin, sackin_check)
  }

})
