#include <vector>
#include <array>
#include <Rcpp.h>

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

// [[Rcpp::export]]
Rcpp::NumericVector Lx_tree_weighted(const Rcpp::IntegerMatrix &edge,
                               const Rcpp::NumericVector &w,
                               const Rcpp::NumericVector &x,
                               int nNodes) {
  Rcpp::NumericVector result(nNodes);

  int nEdges = edge.nrow();
  for (int i = 0; i < nEdges; i++) {

    int u = edge(i, 0) - 1;  // zero-based
    int v = edge(i, 1) - 1;
    double weight = w[i];

    // contribution to Laplacian:
    // L(u) += w * (x[u] - x[v])
    // L(v) += w * (x[v] - x[u])

    double diff_uv = x[u] - x[v];

    result[u] += weight * diff_uv;
    result[v] -= weight * diff_uv;  // = weight * (x[v] - x[u])
  }

  return result;
}
