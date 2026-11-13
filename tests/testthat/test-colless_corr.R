context("colless_corr")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    colless_c <- treestats::colless_corr(extant_tree)
    colless_check <- treebalance::collessI(extant_tree, method = "corrected")
    testthat::expect_equal(colless_c, colless_check)

    # now, using ltable:
    focal_ltab <- treestats::phylo_to_l(extant_tree)

    colless_cl <- treestats::colless_corr(focal_ltab)
    testthat::expect_equal(colless_cl, colless_check)

    ## with extinct lineages:
    colless_c <- treestats::colless_corr(complete_tree)
    colless_check <- treebalance::collessI(complete_tree, method = "corrected")
    testthat::expect_equal(colless_c, colless_check)

    focal_ltab <- treestats::phylo_to_l(complete_tree)

    colless_cl <- treestats::colless_corr(focal_ltab)
    testthat::expect_equal(colless_cl, colless_check)
  }
})

test_that("yule", {
  testthat::skip_on_cran()
  test_yule(treestats::colless_corr)
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {

    test_func <- function(phy) {
      treebalance::collessI(phy, method = "corrected")
    }

    test_polytomies(treestats::colless_corr, test_func)
  }
})


test_that("wrong_object", {
  test_abuse(treestats::colless_corr)
})
