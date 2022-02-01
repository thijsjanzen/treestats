context("crown_age")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  age1 <- treestats::crown_age(focal_tree)
  age2 <- max(adephylo::distRoot(focal_tree))
  testthat::expect_equal(age1, age2)

  # now, more challenging, with extinct lineages:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.3)[[1]]

  age1 <- treestats::crown_age(focal_tree)
  age2 <- max(adephylo::distRoot(focal_tree))
  testthat::expect_equal(age1, age2)

  # and now with root:
  focal_tree <- TreeSim::sim.bd.age(age = 3,
                                    numbsim = 1,
                                    lambda = 1, mu = 0.3,
                                    mrca = TRUE)[[1]]

  age1 <- treestats::crown_age(focal_tree)
  age2 <- max(adephylo::distRoot(focal_tree))
  testthat::expect_equal(age1, age2)

  # and now with root:
  focal_tree <- TreeSim::sim.bd.age(age = 3,
                                    numbsim = 1,
                                    lambda = 1, mu = 0.3,
                                    mrca = FALSE)[[1]]

  age1 <- treestats::crown_age(focal_tree)
  age2 <- max(adephylo::distRoot(focal_tree))
  testthat::expect_equal(age1, age2)

})
