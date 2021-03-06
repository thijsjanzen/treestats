context("gamma")



test_that("usage", {
  if (requireNamespace("castor") &&
      requireNamespace("phytools") &&
      requireNamespace("DDD")) {
    check_gamma <- function(focal_tree) {
      gammast <- treestats::gamma_statistic(focal_tree)
      gammast_check <- ape::gammaStat(focal_tree)
      gammast_check2 <- castor::gamma_statistic(focal_tree)
      gammast_check3 <- phytools::ltt(focal_tree, plot = FALSE,
                                      gamma = TRUE)$gamma

      ltab <- treestats::phylo_to_l(focal_tree)
      gammast_check4 <- treestats::gamma_statistic(ltab)

      testthat::expect_equal(gammast, gammast_check, tolerance = 1e-4)
      testthat::expect_equal(gammast, gammast_check2, tolerance = 1e-4)
      testthat::expect_equal(gammast, gammast_check3, tolerance = 1e-4)
      testthat::expect_equal(gammast, gammast_check4, tolerance = 1e-4)
    }

    set.seed(42)
    bd_tree <- ape::rphylo(n = 101, birth = 1, death = 0)

    check_gamma(bd_tree)

    # DDD tree has all speciation events at the beginning
    # expect strongly negative gamma.
    dd_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 20)$tas
    check_gamma(dd_tree)
    testthat::expect_lt(treestats::gamma_statistic(dd_tree), 0)
  }
})

test_that("abuse", {
  wrong_object <- list("a" = 5)
  testthat::expect_error(
    treestats::colless(wrong_object),
    "input object has to be phylo or ltable")
})
