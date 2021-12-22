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

  std::vector<size_t> interns(Nnode);

  std::vector<double> edge_length = phy["edge.length"];
  Rcpp::NumericMatrix edge = phy["edge"];

  size_t cnt = 0;
  for (size_t i = 0; i < edge_length.size(); ++i) {
    if (edge(i, 1) > n) {
      interns[cnt] = i;
      cnt++;
    }
  }

  std::vector<double> xx(Nnode);

  for (const auto& i : interns) {
    xx[ edge(i, 1) - n - 1 ] = xx[edge(i, 0) - n - 1] + edge_length[i];
  }

  int N = edge_length.size() - 1;
  double depth = xx[edge(N, 0) - n - 1] +  edge_length[N];
  for (auto& i : xx) {
    i = depth - i;
  }
  return xx;
}

#endif
