context("treeness")

test_that("usage", {
  set.seed(42)
  # DDD tree has all speciation events at the beginning
  # expect very low stemminess
  dd_tree1 <- DDD::dd_sim(pars = c(1, 0, 10), age = 20)$tas
  treeness_1 <- treestats::treeness(dd_tree1)
  testthat::expect_lt(treeness_1, 0.5)

  dd_tree1_ltab <- treestats::phylo_to_l(dd_tree1)
  treeness1_ltab <- treestats::treeness(dd_tree1_ltab)
  testthat::expect_equal(treeness_1, treeness1_ltab)


  dd_tree2 <- DDD::dd_sim(pars = c(0.3, 0, 100), age = 20)$tas
  treeness_2 <- treestats::treeness(dd_tree2)
  testthat::expect_gt(treeness_2, treeness_1)

  dd_tree2_ltab <- treestats::phylo_to_l(dd_tree2)
  treeness_2_ltab <- treestats::treeness(dd_tree2_ltab)
  testthat::expect_equal(treeness_2, treeness_2_ltab)
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
