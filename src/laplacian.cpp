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

// [[Rcpp::export]]
Rcpp::NumericMatrix prep_lapl_spec(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);

  auto num_nodes = phy["Nnode"];
  Rcpp::StringVector tips = phy["tip.label"];
  auto num_tips = tips.size();

  std::vector< std::vector< double >> lapl_mat = dist_nodes(edge,
                                                            el,
                                                            num_tips,
                                                            num_nodes);
  Rcpp::NumericMatrix res(lapl_mat.size(), lapl_mat[0].size());

  for (size_t i = 0; i < lapl_mat.size(); ++i) {
    for (size_t j = 0; j < lapl_mat[i].size(); ++j) {
      res(i, j) = lapl_mat[i][j];
    }
    res(i, i) = - std::accumulate(lapl_mat[i].begin(), lapl_mat[i].end(), 0.0);
  }

  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix prep_adj_mat(const std::vector<int>& parent_list,
                                 const std::vector<double>& br_len,
                                 bool use_br_len) {
  auto max_num = *std::max_element(parent_list.begin(), parent_list.end());
  Rcpp::NumericMatrix out_mat(max_num, max_num);

  if (use_br_len) {
    for (size_t i = 0; i < parent_list.size(); i += 2) {
      auto x = parent_list[i] - 1;
      auto y = parent_list[i + 1] - 1;
      out_mat(x, y) = out_mat(y, x) = br_len[i / 2];
    }
  } else {
    for (size_t i = 0; i < parent_list.size(); i += 2) {
      auto x = parent_list[i] - 1;
      auto y = parent_list[i + 1] - 1;
      out_mat(x, y) = out_mat(y, x) = 1;
    }
  }
  return(out_mat);
}


