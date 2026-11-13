context("colless_quad")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    colless_c <- treestats::colless_quad(extant_tree)
    colless_check <- treebalance::collessI(extant_tree, method = "quadratic")
    testthat::expect_equal(colless_c, colless_check)

    # now, using ltable:
    focal_ltab <- treestats::phylo_to_l(extant_tree)

    colless_cl <- treestats::colless_quad(focal_ltab)
    testthat::expect_equal(colless_cl, colless_check)

    ## with extinct lineages:

    colless_c <- treestats::colless_quad(complete_tree)
    colless_check <- treebalance::collessI(complete_tree, method = "quadratic")
    testthat::expect_equal(colless_c, colless_check)

    focal_ltab <- treestats::phylo_to_l(complete_tree)

    colless_cl <- treestats::colless_quad(focal_ltab)
    testthat::expect_equal(colless_cl, colless_check)
  }
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {

    test_func <- function(phy) {
      treebalance::collessI(phy, method = "quadratic")
    }

    test_polytomies(treestats::colless_quad, test_func)
  }
})

test_that("yule", {
  testthat::skip_on_cran()
  test_yule(treestats::colless_quad)
})

test_that("basis equality", {
 # in Bauman & Rosenberg 2026, they show that:
 # QC = 2 * tot_coph + sackin - n * (n - 1), where n is the number of tips

 phy   <- ape::rphylo(n = 100, birth = 1, death = 0)
 res   <- treestats::colless_quad(phy)
 n     <- treestats::number_of_lineages(phy) # 100
 sack  <- treestats::sackin(phy)
 tcoph <- treestats::tot_coph(phy)

 testthat::expect_equal(res,
                        2 * tcoph + sack  - n * (n - 1))
})


test_that("wrong_object", {
  test_abuse(treestats::colless_quad)
})
