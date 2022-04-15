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

#endif
