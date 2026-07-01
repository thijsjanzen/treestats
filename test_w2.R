calcwiener <- function (tree, norm = FALSE, weight = FALSE)
{
  if (is.null(tree$subtreeSizes)) {
    tree = treeCentrality:::addSubtreeSizes(tree)
  }
  q = rowSums(tree$subtreeSizes) + 1
  cat(q, "\n")
  n = length(q)
  N = 2 * n + 1
  stopifnot(q[1] == N)
  if (weight) {
    edges = tree$edge
    lengths = tree$edge.length
    stopifnot(all(lengths >= 0))
    W = 0
    for (ind in 1:nrow(tree$edge)) {
      curEndpoint = edges[ind, 2]
      curQ = ifelse(curEndpoint > (n + 2), q[curEndpoint -
                                               (n + 1)], 1)
      W = W + curQ * (N - curQ) * lengths[ind]
    }
  }
  else {
    W = sum(q * (N - q)) + (N - 1) * (n + 1)
  }
  W = W/ifelse(norm, choose(N, 2), 1)
  W
}

set.seed(5)
extant_tree <- ape::rphylo(n = 5, birth = 1, death = 0)
a1 <- calcwiener(extant_tree, weight = TRUE)
a2 <- treestats::wiener(extant_tree, weight = TRUE)
cat(a1, a2, "\n")