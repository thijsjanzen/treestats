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
#include <numeric>

using ltable = std::vector< std::array<double, 4>>;

double calc_mntd_ltable(const ltable& ltable_) {
  std::vector<double> dist(ltable_.size() + 1, -1);

  for (const auto& i : ltable_) {
    auto parent = std::abs(i[1]);
    auto daughter = std::abs(i[2]);
    auto dist_to_nearest_taxon = 2 * i[0];
    if (i[3] != -1) {
      dist_to_nearest_taxon = i[0] + (i[0] - i[3]);
    }

    if (daughter < 0 || daughter > dist.size()) {
      throw std::out_of_range("daughter outside dist");
    }

    dist[daughter] = dist_to_nearest_taxon;

    if (parent < 0 || parent > dist.size()) {
      throw std::out_of_range("parent outside dist");
    }

    if (dist[parent] > 0) {
      if (dist_to_nearest_taxon < dist[parent]) {
        dist[parent] = dist_to_nearest_taxon;
      }
    } else {
      dist[parent] = dist_to_nearest_taxon;
    }
  }

  dist[0]= 0.0;
  auto sum_dist = std::accumulate(dist.begin(), dist.end(), 0.0);
  return sum_dist * 1.0 / ltable_.size();
}

double calc_mntd_stat(const std::vector< std::array< size_t, 2 >>& edge,
                      const std::vector<double>& el) {
  size_t root_no = edge[0][0];
  size_t max_num = 0;
  for (const auto& i : edge) {
    if (i[0] > max_num) max_num = i[0];
    if (i[0] < root_no) root_no = i[0];
  }

  std::vector<double> node_heights(max_num + 1, 0);
  for (size_t i = 0; i < edge.size(); ++i) {
    node_heights[edge[i][1]] = node_heights[edge[i][0]] + el[i];
  }

  // first N entries are distance from root to tips.
  double crown_age = *std::max_element(node_heights.begin(),
                                       node_heights.begin() + root_no);

  for (auto& i : node_heights) {
    i = crown_age - i;
  }

  // and now we calculate mntd, this is always the distance
  // to the parent root * 2.
  double mntd = 0.0;
  for (const auto& i : edge) {
    if (i[1] < root_no) {  // we now have a tip
      mntd += node_heights[i[0]] * 2;
    }
  }

  mntd *= 1.0 / (root_no - 1);
  return(mntd);
}

// this improved version was fully cooked by chatGPT
double calc_var_mpd_stat(
    const std::vector<std::array<size_t, 2>>& edge,
    const std::vector<double>& el) {

  const size_t n_edges = edge.size();
  const size_t n_nodes = n_edges + 1;
  const size_t n_tips = (n_edges + 2) / 2;

  // Number of tips below each node.
  std::vector<size_t> n(n_nodes, 0);

  // Sum of distances from descendant tips to node.
  std::vector<double> s(n_nodes, 0.0);

  // Sum of squared distances from descendant tips to node.
  std::vector<double> q(n_nodes, 0.0);

  // Sum of pairwise distances within subtree.
  double total_sum = 0.0;

  // Sum of squared pairwise distances within subtree.
  double total_sum_sq = 0.0;

  /*
   * ape's edge matrix is normally ordered such that children occur
   * after their parents. Therefore process edges backwards.
   *
   * First initialize tip counts.
   */
  for (size_t i = 0; i < n_edges; ++i) {
    const size_t child = edge[i][1] - 1;

    // Tips are the first n_tips nodes in ape numbering.
    if (child < n_tips)
      n[child] = 1;
  }

  /*
   * Process each edge from the tips towards the root.
   */
  for (size_t i = n_edges; i-- > 0;) {
    const size_t parent = edge[i][0] - 1;
    const size_t child  = edge[i][1] - 1;
    const double length = el[i];

    const size_t nc = n[child];

    if (nc == 0)
      continue;

    // Distances from child descendants to parent.
    const double sc =
      s[child] + static_cast<double>(nc) * length;

    const double qc =
      q[child]
    + 2.0 * length * s[child]
    + static_cast<double>(nc) * length * length;

    /*
     * Combine this child with all previously processed children
     * of the parent.
     */
    const size_t np = n[parent];

    if (np > 0) {
      total_sum +=
        static_cast<double>(nc) * s[parent]
      + static_cast<double>(np) * sc;

      total_sum_sq +=
      static_cast<double>(nc) * q[parent]
      + static_cast<double>(np) * qc
      + 2.0 * s[parent] * sc;
    }

    n[parent] += nc;
    s[parent] += sc;
    q[parent] += qc;
  }

  const double pairs =
    static_cast<double>(n_tips) *
    static_cast<double>(n_tips - 1) / 2.0;

  if (pairs == 0.0)
    return 0.0;

  // Population variance of all unordered tip pairs.
  const double mean = total_sum / pairs;

  return total_sum_sq / pairs - mean * mean;
}
