context("sackin")

test_that("usage", {
  if (requireNamespace("apTreeshape")) {
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

    testthat::expect_equal(treestats::sackin(focal_tree,
                                             normalization = "yule"),
                           treestats::sackin(ltab,
                                             normalization = "yule"))

    testthat::expect_equal(treestats::sackin(focal_tree, normalization = "pda"),
                           treestats::sackin(ltab, normalization = "pda"))
  }

  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    sackin <- treestats::sackin(focal_tree)
    sackin_check <- treebalance::sackinI(focal_tree)
    testthat::expect_equal(sackin, sackin_check)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::sackin(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::sackin(list()),
    "input object has to be phylo or ltable"
  )
})
