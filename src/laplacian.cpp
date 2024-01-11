// Copyright 2022 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Laplacian matrix statistics
//
//
//
#include <Rcpp.h>
#include "util.h"         // NOLINT [build/include_subdir]
#include "dist_nodes.h"   // NOLINT [build/include_subdir]

//' function to create laplacian matrix
//' @param phy phy
//' @return numericmatrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix prep_lapl_spec(const Rcpp::List& phy) {
  int nnode = phy["Nnode"];
  if (1 + 2 * nnode >= 46340) {   // floor(sqrt(2^31 - 1))
    throw std::runtime_error("tree too big to calculate laplacian spectrum");
  }


  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);

  std::vector< std::vector< double >> lapl_mat = dist_nodes(edge, el);
  Rcpp::NumericMatrix res(lapl_mat.size(), lapl_mat[0].size());

  for (size_t i = 0; i < lapl_mat.size(); ++i) {
    for (size_t j = 0; j < lapl_mat[i].size(); ++j) {
      res(i, j) = lapl_mat[i][j];
    }
    res(i, i) = - std::accumulate(lapl_mat[i].begin(), lapl_mat[i].end(), 0.0);
  }

  return res;
}
