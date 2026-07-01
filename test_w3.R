betw <- function (tree, weight = FALSE)
{
  if (is.null(tree$subtreeSizes)) {
    tree = treeCentrality:::addSubtreeSizes(tree)
  }
  Tab = tree$subtreeSizes
  n = nrow(Tab)
  rSums = rowSums(Tab)
  Centralities = c(rep(0, n + 1), Tab[, 1] * Tab[, 2] + rSums *
                     (2 * n - rSums))
  Centralities
}

test_func <- function(phy) {
  df <- as.data.frame(cbind(phy$edge,
                            weight = phy$edge.length))
  g <- igraph::graph_from_data_frame(df, directed = FALSE)
  ref_betweenness <- igraph::betweenness(g)
  return((ref_betweenness))
}




#extant_tree <- ape::read.tree(text = "((A:0.2,B:0.2):0.5,(C:0.3,D:0.3,E:0.3):0.5);")
# a1 <- betw(extant_tree)
# a2 <- max(a1)
set.seed(5)
extant_tree <- ape::rphylo(n = 5, 1, 0)
plot(extant_tree)
a3 <- treestats::max_betweenness(extant_tree)

a4 <- test_func(extant_tree)
cat("\n")
cat(a4, "\n")
cat("\n")
cat("results: ", a3, max(a4), "\n")

if (1 == 2) {
set.seed(5)
extant_tree <- ape::rphylo(n = 5, birth = 1, death = 0)
a1 <- betw(extant_tree)
cat(a1, "\n")
a2 <- max(a1)
a3 <- treestats::max_betweenness(extant_tree)
cat(a2, a3, "\n")
a4 <- test_func(extant_tree)
a4
}
