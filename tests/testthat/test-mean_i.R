context("mean I")

test_that("usage", {
  if (requireNamespace("treebalance")) {

    test_func <- function(phy) {
      return(treebalance::IbasedI(phy, method = "mean",
                           correction = "prime", logs = FALSE))
    }

    standard_test(treestats::mean_i, test_func)

    test_polytomies(treestats::mean_i, test_func)
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
