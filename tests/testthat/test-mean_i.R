context("mean I")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::mean_i(focal_tree)
    a2 <- treebalance::IbasedI(focal_tree, method = "mean",
                               correction = "prime", logs = FALSE)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::mean_i(focal_tree),
                           treestats::mean_i(ltab))


    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    a1 <- treestats::mean_i(focal_tree)
    a2 <- treebalance::IbasedI(focal_tree, method = "mean",
                               correction = "prime", logs = FALSE)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::mean_i(focal_tree),
                           treestats::mean_i(ltab))

    av <- treestats::i_stat(focal_tree)
    testthat::expect_equal(a1, av)
    testthat::expect_equal(treestats::i_stat(focal_tree),
                           treestats::i_stat(ltab))

  }
})

test_that("abuse", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 3, birth = 1, death = 0)
  testthat::expect_warning(treestats::mean_i(focal_tree),
              "I statistic is only available for trees with at least 4 tips.")

  focal_ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_warning(treestats::mean_i(focal_ltab),
              "I statistic is only available for trees with at least 4 tips.")
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::mean_i(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mean_i(list()),
    "input object has to be phylo or ltable"
  )
})
