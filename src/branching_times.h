#ifndef BRANCHING_TIMES_H
#define BRANCHING_TIMES_H

#include <Rcpp.h>

//' function to calculate branching times of a tree, approx 20 times faster
//' than the ape equivalent.
//' @param phy phylo object
//' @export
// [[Rcpp::export]]
std::vector< double > branching_times(const Rcpp::List& phy) {

  size_t Nnode = phy["Nnode"];
  size_t n = Nnode + 1;

  std::vector<double> edge_length = phy["edge.length"];
  Rcpp::NumericMatrix edge = phy["edge"];

  std::vector<double> xx(Nnode, 0.0);

  for (size_t i = 0; i < edge_length.size(); ++i) {
    if (edge(i, 1) > n) {
      auto target_index = edge(i, 1) - n - 1; // -1 because of R to C++ indexing
      auto source_index = edge(i, 0) - n - 1;

      if (target_index < 0 || target_index >= xx.size()) { // range checks
        throw std::out_of_range("target_index out of bounds");
      }

      if (source_index < 0 || source_index >= xx.size()) {
        throw std::out_of_range("source_index out of bounds");
      }

      xx[ target_index  ] = xx[ source_index] + edge_length[i];
    }
  }

  auto edge_index = edge(edge_length.size() - 1, 0) - n - 1;
  if (edge_index < 0 || edge_index >= xx.size()) {
    throw std::out_of_range("edge_index out of bounds");
  }

  double depth = xx[edge_index] +  edge_length[edge_length.size() - 1];
  for (auto& i : xx) {
    i = depth - i;
  }
  return xx;
}

#endif
