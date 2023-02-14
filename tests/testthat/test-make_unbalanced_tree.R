context("make_unbalanced_tree")

test_that("usage", {
  set.seed(42)
  tree_size <- 8
  focal_tree <- ape::rphylo(n = tree_size, birth = 1, death = 0)

  num_steps <- tree_size * 2
  for (method in c("youngest", "random",
                   "terminal-random", "terminal-youngest")) {
    unbal_tree <- treestats::make_unbalanced_tree(focal_tree,
                                                  unbal_steps = num_steps,
                                                  method = method)
    max_d <- treestats::max_depth(unbal_tree)
    testthat::expect_equal(max_d, tree_size - 1)
  }
})
