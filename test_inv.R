set.seed(42)
focal_tree <- ape::rphylo(n = 3, 1, 0)
focal_tree <- stats::reorder(focal_tree)
focal_tree$edge.length <- c(0.5, 1, 1, 1.5)
expected_vals <- c(1 + 1 / 0.5, 1 + 1 / 0.5, 1 / 1.5)

inv_dist <- treestats::inv_branch_dist(focal_tree, add_one = FALSE)

cat("treestats: ", inv_dist, "\n")
cat("expected: ", expected_vals, "\n")

testthat::expect_equal(as.vector(inv_dist), expected_vals)

ltab <- treestats::phylo_to_l(focal_tree)
inv_dist_ltab <- treestats::inv_branch_dist(ltab, add_one = FALSE)
testthat::expect_equal(sum(as.vector(inv_dist_ltab)), sum(expected_vals))

# now we add one
expected_vals <- c(1 / (1 + 1) + 1 / 1.5, 1 / (1 + 1) + 1 / 1.5, 1 / 2.5)

inv_dist <- treestats::inv_branch_dist(focal_tree, add_one = TRUE)
cat(inv_dist, "\n")
cat(expected_vals, "\n")
testthat::expect_equal(as.vector(inv_dist), expected_vals)

# for (focal_tree in poly_trees) {
#  local_stat <- try_stat(focal_tree, treestats::inv_branch_dist)
#  testthat::expect_true(is.na(local_stat))
# }

