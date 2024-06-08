context("tot_path")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::tot_path_length(focal_tree)
    a2 <- treebalance::totPathLen(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::tot_path_length(focal_tree),
                           treestats::tot_path_length(ltab))


    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    a1 <- treestats::tot_path_length(focal_tree)
    a2 <- treebalance::totPathLen(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::tot_path_length(focal_tree),
                           treestats::tot_path_length(ltab))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::tot_path_length(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::tot_path_length(list()),
    "input object has to be phylo or ltable"
  )
})
