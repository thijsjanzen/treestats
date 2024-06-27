context("pigot_rho")

test_that("usage", {
  set.seed(42)

  if (requireNamespace("DDD")) {
    # DDD tree expected to slow down diversification
    focal_tree <- DDD::dd_sim(pars = c(1, 0, 10), age = 7)$tes
    rho <- treestats::pigot_rho(focal_tree)
    testthat::expect_lt(rho, 0)
  }

  set.seed(42)
  # BD tree, expected increase in diversification
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
  rho <- treestats::pigot_rho(focal_tree)
  testthat::expect_gt(rho, 0)

  ltab <- treestats::phylo_to_l(focal_tree)
  rho2 <- treestats::pigot_rho(ltab)
  testthat::expect_equal(rho, rho2)

  if (requireNamespace("geiger")) {
    # and we can do with extinct trees as well
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2,
                              fossils = TRUE)
    rho4 <- treestats::pigot_rho(focal_tree)
    extant_tree <- geiger::drop.extinct(focal_tree)
    rho5 <- treestats::pigot_rho(extant_tree)
    testthat::expect_false(rho4 == rho5)
  }

  # use very small tree to trigger use of complete method:
  focal_tree <- ape::rphylo(n = 6, birth = 1, death = 0)
  rho <- treestats::pigot_rho(focal_tree)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::pigot_rho(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::pigot_rho(list()),
    "input object has to be phylo or ltable"
  )
})
