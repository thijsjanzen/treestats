context("treeness")

test_that("usage", {
  set.seed(42)
  # DDD tree has all speciation events at the beginning
  # expect very low stemminess
  dd_tree1 <- DDD::dd_sim(pars = c(1, 0, 10), age = 20)$tas
  vv1 <- treestats::treeness(dd_tree1)
  testthat::expect_lt(vv, 0.5)

  dd_tree2 <- DDD::dd_sim(pars = c(0.3, 0, 100), age = 20)$tas
  vv2 <- treestats::treeness(dd_tree2)
  testthat::expect_gt(vv2, vv1)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::treeness(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::treeness(list()),
    "input object has to be phylo or ltable"
  )
})
