#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector phylo_laplacian_matvec_rcpp(
    Rcpp::IntegerMatrix edges,
    Rcpp::NumericVector lengths,
    Rcpp::NumericVector x) {
  // edges: (n_edges × 2) matrix of parent/child node indices (1-based)
  // lengths: vector of branch lengths (same length as n_edges)
  // x: vector of node values
  // output: y = L * x

  int n_nodes = x.size();
  Rcpp::NumericVector y(n_nodes);
  int n_edges = edges.nrow();

  for(int e = 0; e < n_edges; e++) {
    int i = edges(e, 0) - 1; // parent
    int j = edges(e, 1) - 1; // child
    double w = 1.0 / lengths[e]; // weight = 1 / branch length

    // contributions for both directions (undirected)
    y[i] += w * (x[i] - x[j]);
    y[j] += w * (x[j] - x[i]);
  }
  return y;
}

// [[Rcpp::export]]
Rcpp::NumericVector Ax_tree(const Rcpp::IntegerMatrix &edge,
                            const Rcpp::NumericVector &lengths,
                            const Rcpp::NumericVector &x,
                            int nNodes) {
  Rcpp::NumericVector result(nNodes);

  int nEdges = edge.nrow();
  for (int i = 0; i < nEdges; i++) {
    int u = edge(i, 0) - 1; // convert to zero-based
    int v = edge(i, 1) - 1;

    double w = lengths[i];

    // adjacency contribution
    result[u] += w * x[v];
    result[v] += w * x[u];
  }

  return result;
}
