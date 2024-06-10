context("minmax_adj")

test_that("usage", {
  if (requireNamespace("igraph")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 5, birth = 1, death = 0)

    a1_1 <- treestats::minmax_adj(focal_tree)

    df <- as.data.frame(cbind(focal_tree$edge,
                             weight = focal_tree$edge.length))
    g <- igraph::graph_from_data_frame(df, directed = FALSE)

    adj_mat <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
    ref <- eigen(adj_mat)$values

    testthat::expect_equal(a1_1$min, min(ref[ref > 0]))
    testthat::expect_equal(a1_1$max, max(ref))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::minmax_adj(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::minmax_adj(list()),
    "input object has to be phylo or ltable"
  )
})
