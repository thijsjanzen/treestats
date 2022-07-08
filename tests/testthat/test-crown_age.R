context("crown_age")

test_that("usage", {
  if (requireNamespace("adephylo") &&
      requireNamespace("TreeSim")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    age1 <- treestats::tree_height(focal_tree)
    age2 <- max(adephylo::distRoot(focal_tree))
    testthat::expect_equal(age1, age2)

    # now, more challenging, with extinct lineages:
    focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                       numbsim = 1,
                                       lambda = 1, mu = 0.3)[[1]]

    age1 <- treestats::tree_height(focal_tree)
    age2 <- max(adephylo::distRoot(focal_tree))
    testthat::expect_equal(age1, age2)

    # and now with crown age:
    focal_tree <- TreeSim::sim.bd.age(age = 3,
                                      numbsim = 1,
                                      lambda = 1, mu = 0.3,
                                      mrca = TRUE)[[1]]

    age1 <- treestats::crown_age(focal_tree)
    testthat::expect_equal(age1, 3)


    # and now without root.edge
    focal_tree <- ape::rphylo(100, 1, 0)
    age1 <- treestats::crown_age(focal_tree)
    age2 <- max(adephylo::distRoot(focal_tree))
    testthat::expect_equal(age1, age2)
  }
})
