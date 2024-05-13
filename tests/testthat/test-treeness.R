context("treeness")

test_that("usage", {
  set.seed(42)
  # DDD tree has all speciation events at the beginning
  # expect very low stemminess
  dd_tree1 <- DDD::dd_sim(pars = c(1, 0, 10), age = 20)$tas
  vv1 <- treestats::treeness(dd_tree1)
  testthat::expect_lt(vv, 0.5)

  dd_tree1_ltab <- treestats::phylo_to_l(dd_tree1)
  vv1_ltab <- treestats::treeness(dd_tree1_ltab)
  testthat::expect_equal(vv1, vv1_ltab)


  dd_tree2 <- DDD::dd_sim(pars = c(0.3, 0, 100), age = 20)$tas
  vv2 <- treestats::treeness(dd_tree2)
  testthat::expect_gt(vv2, vv1)

  dd_tree2_ltab <- treestats::phylo_to_l(dd_tree2)
  vv2_ltab <- treestats::treeness(dd_tree2_ltab)
  testthat::expect_equal(vv2, vv2_ltab)
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
