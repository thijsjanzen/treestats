#ifndef util_h
#define util_h

#include <vector>
#include <array>
#include <Rcpp.h>

using ltable = std::vector< std::array<double, 4>>;
using edge_table = std::vector< std::array< size_t, 2 >>;

// inline to allow inclusion in multiple cpp
inline ltable convert_to_ltable(const Rcpp::NumericMatrix& mat_in) {
  ltable out(mat_in.nrow());

  for (size_t i = 0; i < mat_in.nrow(); ++i) {
    std::array<double, 4> row_entry = {mat_in(i, 0), mat_in(i, 1),
                                       mat_in(i, 2), mat_in(i, 3) };
    out[i] = row_entry;
  }
  return out;
}


// short util functions:
inline edge_table phy_to_edge(const Rcpp::List& phy) {
  Rcpp::NumericMatrix edge = phy["edge"];
  edge_table local_edge(edge.nrow());
  for (size_t i = 0; i < edge.nrow(); ++i) {
    local_edge[i] = {static_cast<size_t>(edge(i, 0)),
                     static_cast<size_t>(edge(i, 1))};
  }
  return local_edge;
}

inline std::vector<double> phy_to_el(const Rcpp::List& phy) {
  Rcpp::NumericVector el = phy["edge.length"];
  std::vector<double> el_cpp(el.begin(), el.end());
  return el_cpp;
}

#endif
