context("gamma")



test_that("usage", {
  if (requireNamespace("castor") &&
      requireNamespace("DDD")) {
    check_gamma <- function(focal_tree) {
      gammast <- treestats::gamma_statistic(focal_tree)
      gammast_check <- ape::gammaStat(focal_tree)
      gammast_check2 <- castor::gamma_statistic(focal_tree)

      ltab <- treestats::phylo_to_l(focal_tree)
      gammast_check4 <- treestats::gamma_statistic(ltab)

      testthat::expect_equal(gammast, gammast_check, tolerance = 1e-4)
      testthat::expect_equal(gammast, gammast_check2, tolerance = 1e-4)
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

test_that("wrong_object", {
  testthat::expect_error(
    treestats::gamma_statistic(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::gamma_statistic(list()),
    "input object has to be phylo or ltable"
  )
})
