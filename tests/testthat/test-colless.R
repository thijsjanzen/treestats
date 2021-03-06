context("colless")

test_that("usage", {
  if (requireNamespace("apTreeshape")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    colless <- treestats::colless(focal_tree)
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree))
    testthat::expect_equal(colless, colless_check)



    colless <- treestats::colless(focal_tree, normalization = "yule")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                          norm = "yule")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)

    colless <- treestats::colless(focal_tree, normalization = "pda")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                          norm = "pda")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)


    # now, using ltable:
    ffocal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
    focal_ltab <- treestats::phylo_to_l(focal_tree)

    colless <- treestats::colless(focal_ltab)
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree))
    testthat::expect_equal(colless, colless_check)

    colless <- treestats::colless(focal_ltab, normalization = "yule")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                          norm = "yule")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)

    colless <- treestats::colless(focal_ltab, normalization = "pda")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                          norm = "pda")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)


    ## with extinct lineages:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    colless <- treestats::colless(focal_tree)
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree))
    testthat::expect_equal(colless, colless_check)

    colless <- treestats::colless(focal_tree, normalization = "yule")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                          norm = "yule")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)

    colless <- treestats::colless(focal_tree, normalization = "pda")
    colless_check <- apTreeshape::colless(apTreeshape::as.treeshape(focal_tree),
                                          norm = "pda")
    testthat::expect_equal(colless, colless_check, tol = 1e-5)
  }
})

test_that("abuse", {
  wrong_object <- list("a" = 5)
  testthat::expect_error(
    treestats::colless(wrong_object),
    "input object has to be phylo or ltable")
})
