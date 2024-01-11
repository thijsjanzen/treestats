context("L2phylo")

test_that("usage", {
  if (requireNamespace("geiger") &&
      requireNamespace("DDD")) {
    set.seed(42)

    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    ltable_1 <- treestats::phylo_to_l(focal_tree)

    tree2 <- treestats::l_to_phylo(ltable_1)
    tree3 <- DDD::L2phylo(ltable_1)

    ax <- geiger::is.extinct(tree2)
    ay <- geiger::is.extinct(tree3)

    testthat::expect_true(all.equal(ax, ay))

    diff_edge_length <- sort(tree2$edge.length) - sort(tree3$edge.length)
    testthat::expect_equal(mean(diff_edge_length), 0.0, tolerance = 0.001)

    # and now with a tree with extinct tips
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.3, fossils = TRUE)

    ltable_1 <- treestats::phylo_to_l(focal_tree)

    tree2 <- treestats::l_to_phylo(ltable_1)
    tree3 <- DDD::L2phylo(ltable_1)

    ax <- geiger::is.extinct(tree2)
    ay <- geiger::is.extinct(tree3)

    testthat::expect_true(all.equal(ax, ay))
    diff_edge_length <- sort(tree2$edge.length) - sort(tree3$edge.length)
    testthat::expect_equal(mean(diff_edge_length), 0.0, tolerance = 0.001)

    # and again, but retain extinct lineages
    tree2 <- treestats::l_to_phylo(ltable_1, drop_extinct = FALSE)
    tree3 <- DDD::L2phylo(ltable_1, dropextinct = FALSE)

    ax <- geiger::is.extinct(tree2)
    ay <- geiger::is.extinct(tree3)

    testthat::expect_true(all.equal(ax, ay))
    diff_edge_length <- sort(tree2$edge.length) - sort(tree3$edge.length)
    testthat::expect_equal(mean(diff_edge_length), 0.0, tolerance = 0.001)
  }
})

test_that("newick", {
  if (requireNamespace("geiger")) {
    set.seed(42)
    focal_tree <- ape::rbdtree(birth = 1, death = 0, Tmax = 5)

    ltable_1 <- treestats::phylo_to_l(focal_tree)

    newick_str <- treestats::ltable_to_newick(ltable_1)
    tree_2 <- ape::read.tree(text = newick_str)
    v1 <- as.vector(unlist(treestats::calc_all_stats(focal_tree)))
    v2 <- as.vector(unlist(treestats::calc_all_stats(tree_2)))
    testthat::expect_equal(v1, v2, tolerance = 1e-3)

    focal_tree <- ape::rphylo(n = 100,
                              birth = 1, death = 0.2, T0 = 2, fossils = TRUE)

    ltable_1 <- treestats::phylo_to_l(focal_tree)

    newick_str <- treestats::ltable_to_newick(ltable_1)
    tree_2 <- ape::read.tree(text = newick_str)
    tree_2 <- geiger::drop.extinct(tree_2)
    focal_tree <- geiger::drop.extinct(focal_tree)
    v1 <- as.vector(unlist(treestats::calc_all_stats(focal_tree)))
    v2 <- as.vector(unlist(treestats::calc_all_stats(tree_2)))
    testthat::expect_equal(v1, v2, tolerance = 1e-3)
  }
})
