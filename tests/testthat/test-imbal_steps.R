context("imbal_steps")

test_that("usage", {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 32, birth = 1, death = 0)
    brts <- treestats::branching_times(focal_tree)
    if (requireNamespace("nodeSub")) {
      focal_tree <- nodeSub::create_balanced_tree(brts)
      a1 <- treestats::imbalance_steps(focal_tree)
      a2 <- 32 - log(32, 2) - 1
      testthat::expect_equal(a1, a2)

      ltab <- treestats::phylo_to_l(focal_tree)
      testthat::expect_equal(treestats::imbalance_steps(focal_tree),
                             treestats::imbalance_steps(ltab))

      # normalization
      a1 <- treestats::imbalance_steps(focal_tree,
                                       normalization = TRUE)
      testthat::expect_equal(a1, 1.0)
    }
})

test_that("abuse", {
  tree1 <- ape::rphylo(n = 2, birth = 1, death = 0)
  testthat::expect_message(
    treestats::imbalance_steps(tree1, normalization = TRUE)
  )

  testthat::expect_message(
    treestats::imbalance_steps(treestats::phylo_to_l(tree1),
                               normalization = TRUE)
  )
})


test_that("wrong_object", {
  testthat::expect_error(
    treestats::imbalance_steps(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::imbalance_steps(list()),
    "input object has to be phylo or ltable"
  )
})
