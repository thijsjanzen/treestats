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
#pragma once

#include <vector>
#include <array>
#include <Rcpp.h>

using ltable = std::vector< std::array<double, 4>>;
using edge_table = std::vector< std::array< size_t, 2 >>;

// inline to allow inclusion in multiple cpp
inline ltable convert_to_ltable(const Rcpp::NumericMatrix& mat_in) {
  ltable out(mat_in.nrow());

  for (int i = 0; i < mat_in.nrow(); ++i) {
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
  for (int i = 0; i < edge.nrow(); ++i) {
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
