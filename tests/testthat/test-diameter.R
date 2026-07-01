context("diameter")

test_that("usage", {
  a1_1 <- treestats::diameter(extant_tree, weight = TRUE)
  a1_2 <- treestats::diameter(extant_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 8.750583  # treeCentrality::computeDiameter(focal_tree,
  #                                                   weight = TRUE))
  a2_2 <- 21        # treeCentrality::computeDiameter(focal_tree,
  #                                                      weight = FALSE))

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)

  ltab <- treestats::phylo_to_l(extant_tree)
  testthat::expect_equal(treestats::diameter(extant_tree, weight = TRUE),
                         treestats::diameter(ltab, weight = TRUE))

  testthat::expect_equal(treestats::diameter(extant_tree, weight = FALSE),
                         treestats::diameter(ltab, weight = FALSE))

  a1_1 <- treestats::diameter(complete_tree, weight = TRUE)
  a1_2 <- treestats::diameter(complete_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 11.21469 # treeCentrality::computeDiameter(focal_tree,
  #                                                  weight = TRUE))
  a2_2 <- 21       # treeCentrality::computeDiameter(focal_tree,
  #                                                  weight = FALSE))


  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)
  testthat::expect_equal(a1_2, a2_2)

  ltab <- treestats::phylo_to_l(complete_tree)
  testthat::expect_equal(treestats::diameter(complete_tree, weight = TRUE),
                         treestats::diameter(ltab,       weight = TRUE))

  testthat::expect_equal(treestats::diameter(complete_tree, weight = FALSE),
                         treestats::diameter(ltab,       weight = FALSE))
})

test_that("polytomies", {
  if (requireNamespace("igraph")) {

    test_func <- function(phy) {
      df <- as.data.frame(cbind(phy$edge,
                                weight = phy$edge.length))
      g <- igraph::graph_from_data_frame(df, directed = FALSE)
      return(igraph::diameter(g))
    }
    # TODO: fix this code!
    # test_polytomies(treestats::diameter, test_func)
  }
})


test_that("wrong_object", {
  testthat::expect_error(
    treestats::diameter(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::diameter(list()),
    "input object has to be phylo or ltable"
  )
})
