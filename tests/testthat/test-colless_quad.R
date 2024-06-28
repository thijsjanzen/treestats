context("colless_quad")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    colless_c <- treestats::colless_quad(focal_tree)
    colless_check <- treebalance::collessI(focal_tree, method = "quadratic")
    testthat::expect_equal(colless_c, colless_check)

    # now, using ltable:
    focal_ltab <- treestats::phylo_to_l(focal_tree)

    colless_cl <- treestats::colless_quad(focal_ltab)
    testthat::expect_equal(colless_cl, colless_check)

    ## with extinct lineages:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    colless_c <- treestats::colless_quad(focal_tree)
    colless_check <- treebalance::collessI(focal_tree, method = "quadratic")
    testthat::expect_equal(colless_c, colless_check)

    focal_ltab <- treestats::phylo_to_l(focal_tree)

    colless_cl <- treestats::colless_quad(focal_ltab)
    testthat::expect_equal(colless_cl, colless_check)
  }
})

test_that("yule", {
  testthat::skip_on_cran()
  found <- c()
  found2 <- c()
  for (r in 1:100) {
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
    found[r] <- treestats::colless_quad(focal_tree, "yule")
    found2[r] <- treestats::colless_quad(treestats::phylo_to_l(focal_tree),
                                         "yule")
  }
  testthat::expect_equal(mean(found), 1.0, tolerance = 0.1)
  testthat::expect_equal(mean(found), mean(found2))
})


test_that("wrong_object", {
  testthat::expect_error(
    treestats::colless_quad(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::colless_quad(list()),
    "input object has to be phylo or ltable"
  )
})
