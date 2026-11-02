context("number_of_lineages")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  num_lin <- treestats::number_of_lineages(focal_tree)

  testthat::expect_equal(num_lin, 100)

  focal_ltab <- treestats::phylo_to_l(focal_tree)
  num_lin2 <- treestats::number_of_lineages(focal_ltab)
  testthat::expect_equal(num_lin, 100)

  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2,
                            fossils = TRUE)

  num_lin <- treestats::number_of_lineages(focal_tree)

  testthat::expect_gte(num_lin, 100)

  focal_ltab <- treestats::phylo_to_l(focal_tree)
  num_lin2 <- treestats::number_of_lineages(focal_ltab)
  testthat::expect_gte(num_lin, 100)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::number_of_lineages(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::number_of_lineages(list()),
    "input object has to be phylo or ltable"
  )
})
