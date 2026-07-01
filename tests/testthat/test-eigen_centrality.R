context("eigenvector")

test_that("usage", {
  a1_1 <- treestats::eigen_centrality(extant_tree, weight = TRUE, scale = FALSE)
  a1_2 <- treestats::eigen_centrality(extant_tree, weight = FALSE, scale = FALSE)



  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 3.585569  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = TRUE))
  a2_2 <- 2.609486  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = FALSE))

  testthat::expect_equal(a1_1$eigenvalue, a2_1, tolerance = 1e-4)
  testthat::expect_equal(a1_2$eigenvalue, a2_2, tolerance = 1e-4)
  a1_3 <- treestats::eigen_centrality(extant_tree, weight = TRUE, scale = TRUE)
  testthat::expect_equal(a1_1$eigenvalue,
                         a1_3$eigenvalue)
  testthat::expect_equal(max(a1_3$eigenvector), 1)

  ltab <- treestats::phylo_to_l(extant_tree)
  testthat::expect_equal(
    treestats::eigen_centrality(extant_tree, weight = TRUE)$eigenvalue,
    treestats::eigen_centrality(ltab, weight = TRUE)$eigenvalue)

  a1_1 <- treestats::eigen_centrality(complete_tree, weight = TRUE)
  a1_2 <- treestats::eigen_centrality(complete_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 3.718206  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = TRUE))
  a2_2 <- 2.631605  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = FALSE))

  testthat::expect_equal(a1_1$eigenvalue, a2_1, tolerance = 1e-4)
  testthat::expect_equal(a1_2$eigenvalue, a2_2, tolerance = 1e-4)

  if (requireNamespace("RSpectra")) {
    a1_1 <- treestats::eigen_centrality(complete_tree, weight = TRUE,
                                        use_rspectra = TRUE)
    a1_2 <- treestats::eigen_centrality(complete_tree, weight = FALSE,
                                        use_rspectra = TRUE)
    testthat::expect_equal(a1_1$eigenvalue, a2_1, tolerance = 1e-4)
    testthat::expect_equal(a1_2$eigenvalue, a2_2, tolerance = 1e-4)
  }

  ltab <- treestats::phylo_to_l(complete_tree)
  testthat::expect_equal(
    treestats::eigen_centrality(complete_tree, weight = TRUE)$eigenvalue,
    treestats::eigen_centrality(ltab,          weight = TRUE)$eigenvalue)

  if (requireNamespace("igraph")) {
    df <- as.data.frame(cbind(extant_tree$edge,
                              weight = extant_tree$edge.length))
    g <- igraph::graph_from_data_frame(df, directed = FALSE)

    ref <- igraph::eigen_centrality(g)

    a1 <- treestats::eigen_centrality(extant_tree)

    testthat::expect_equal(a1$eigenvalue, ref$value)

    ltab <- treestats::phylo_to_l(extant_tree)
    a1_2 <- treestats::eigen_centrality(ltab)
    testthat::expect_equal(a1_2$eigenvalue, ref$value, tolerance = 0.01)
  }

  # compare namespaces
  if (requireNamespace("Matrix")) {
    a1_1 <- treestats::eigen_centrality(extant_tree,
                                        weight = TRUE,
                                        scale = FALSE)
    a2_1 <- treestats::eigen_centrality(extant_tree,
                                        weight = FALSE,
                                        scale = FALSE)

    testthat::with_mocked_bindings(
      {
        # Now `myfun()` should behave as if `data.tree` is not installed
        a1_2 <- treestats::eigen_centrality(extant_tree,
                                            weight = TRUE,
                                            scale = FALSE)
        testthat::expect_equal(a1_1, a1_2)

        a1_3 <- treestats::eigen_centrality(extant_tree,
                                            weight = FALSE,
                                            scale = FALSE)
        testthat::expect_equal(a2_1, a1_3)
      },
      requireNamespace = function(pkg, quietly = TRUE) {
        if (pkg == "Matrix") {
          return(FALSE)
        }
        if (pkg == "RSpectra") {
          return(FALSE)
        }
        # Call the real `requireNamespace` for other packages
        base::requireNamespace(pkg, quietly = TRUE)
      },
      .package = "base"
    )
  }
})

test_that("polytomies", {
  if (requireNamespace("igraph")) {

    test_func <- function(phy) {
      df <- as.data.frame(cbind(phy$edge,
                                weight = phy$edge.length))
      g <- igraph::graph_from_data_frame(df, directed = FALSE)
      ref <- igraph::eigen_centrality(g)
      return(ref$value)
    }

    treestats_func <- function(phy) {
      res <- treestats::eigen_centrality(phy)
      return(res$eigenvalue)
    }

    test_polytomies(treestats_func, test_func)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::eigen_centrality(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::eigen_centrality(list()),
    "input object has to be phylo or ltable"
  )
})
