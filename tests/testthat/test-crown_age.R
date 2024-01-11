context("crown_age")

test_that("usage", {
  if (requireNamespace("adephylo")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    age1 <- treestats::tree_height(focal_tree)
    age2 <- max(adephylo::distRoot(focal_tree))
    testthat::expect_equal(age1, age2)

    # now, more challenging, with extinct lineages:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.3, fossils = TRUE)

    age1 <- treestats::tree_height(focal_tree)
    age2 <- max(adephylo::distRoot(focal_tree))
    testthat::expect_equal(age1, age2)

    # and now with crown age:
    focal_tree <- ape::rbdtree(birth = 1, death = 0.3, Tmax = 3)

    age1 <- treestats::crown_age(focal_tree)
    testthat::expect_equal(age1, 3)

    # and now without root.edge
    focal_tree <- ape::rphylo(100, 1, 0)
    age1 <- treestats::crown_age(focal_tree)
    age2 <- max(adephylo::distRoot(focal_tree))
    testthat::expect_equal(age1, age2)

    age3 <- treestats::crown_age(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(age3, age1)
  }
})


test_that("wrong_object", {
  testthat::expect_error(
    treestats::crown_age(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::crown_age(list()),
    "input object has to be phylo or ltable"
  )
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::tree_height(10),
    "input object has to be of class phylo"
  )

  testthat::expect_error(
    treestats::tree_height(list()),
    "input object has to be of class phylo"
  )
})
