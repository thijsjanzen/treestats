// Copyright 2022 - 2025 Thijs Janzen
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
#include <algorithm>
#include <RcppArmadillo.h>

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

struct entry {
  std::array<size_t, 2> ed;
  double bl;
};

inline void add_entry(std::vector<entry>& cladewise_result,
                  std::vector<bool>& added,
                  const std::vector<entry>& unsorted,
                  size_t id,
                  size_t start) {
  for (size_t i = start; i < unsorted.size(); ++i) {
    auto focal_id = unsorted[i].ed[0];
    if (focal_id == id && added[i] == false) {
      cladewise_result.push_back(unsorted[i]);
      added[i] = true;
      auto daughter = unsorted[i].ed[1];
      add_entry(cladewise_result, added, unsorted, daughter, i);
    }
  }
}


inline void sort_edge_and_edgelength(std::vector< std::array<size_t, 2 >>* edge,
                                     std::vector<double>* edge_length) {
  if ((*edge).size() != (*edge_length).size()) {
    throw std::runtime_error("size mismatch");
  }

  std::vector<entry> everything((*edge).size());
  for (size_t i = 0; i < (*edge).size(); ++i) {
    everything[i].bl = (*edge_length)[i];
    everything[i].ed = (*edge)[i];
  }

  std::sort(everything.begin(), everything.end(),
            [&](auto a, auto b)
            {return a.ed[0] < b.ed[0];});

  // now place back
  for (size_t i = 0; i < everything.size(); ++i) {
    (*edge)[i] = everything[i].ed;
    (*edge_length)[i] = everything[i].bl;
  }
  return;
}
