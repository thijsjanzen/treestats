context("make_unbalanced_tree")

test_that("usage", {
  set.seed(42)
  tree_size <- 8
  focal_tree <- ape::rphylo(n = tree_size, birth = 1, death = 0)

  num_steps <- tree_size * 2
  for (g_method in c("any", "terminal")) {
    for (s_method in c("oldest", "youngest", "random")) {
      unbal_tree <- treestats::make_unbalanced_tree(init_tree = focal_tree,
                                                    unbal_steps = num_steps,
                                                    group_method = g_method,
                                                    selection_method = s_method)
      max_d <- treestats::max_depth(unbal_tree)
      testthat::expect_equal(max_d, tree_size - 1)
    }
  }
})

test_that("abuse", {
  set.seed(42)
  tree_size <- 8
  focal_tree <- ape::rphylo(n = tree_size, birth = 1, death = 0)

  testthat::expect_error(
  unbal_tree <- treestats::make_unbalanced_tree(init_tree = focal_tree,
                                                unbal_steps = num_steps,
                                                group_method = "x",
                                                selection_method = "random"),
  "group method unknown, pick from 'any' and 'terminal'")

  testthat::expect_error(
    unbal_tree <- treestats::make_unbalanced_tree(init_tree = focal_tree,
                                                  unbal_steps = num_steps,
                                                  group_method = "any",
                                                  selection_method = "x"),
    "selection method unknown, pick from 'youngest', 'oldest' or 'random"
  )
})
