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

#include <algorithm>   // std::sort
#include <vector>
#include "phylotree.h"  // NOLINT [build/include_subdir]

namespace width {

class width_tree {
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    int depth;

    node_t() {
      depth = 0;
    }

    void set_depth(size_t parent_depth) {
      depth = 1 + parent_depth;

      if (daughter1 && daughter2) {
        daughter1->set_depth(depth);
        daughter2->set_depth(depth);
      }

      return;
    }
  };

  phylo_tree_t<node_t> tree;
  int root_no;

 public:
  explicit width_tree(const std::vector< int >& tree_edge)
    : tree(make_phylo_tree<node_t, true>(tree_edge)) {
    root_no = tree_edge[0];
    for (size_t i = 2; i < tree_edge.size(); i+=2) {
      if (tree_edge[i] < root_no) root_no = tree_edge[i];
    }

    tree[root_no].set_depth(-1);
  }

  int calc_max_width() {
    std::vector<int> depths(tree.size(), 0);
    for (auto i = tree.begin() + 1; i < tree.end(); ++i) {
      depths[ (*i).depth ]++;
    }
    return *std::max_element(depths.begin(), depths.end());
  }

  int calc_max_del_width() {
    std::vector<int> depths(tree.size(), 0);
    for (auto i = tree.begin() + 1; i < tree.end(); ++i) {
      depths[ (*i).depth ]++;
    }
    std::vector<int> dW(depths.size() - 1);
    for (size_t i = 1; i < depths.size(); ++i) {
      dW[i - 1] = depths[i] - depths[i - 1];
    }
    return(*std::max_element(dW.begin(), dW.end()));
  }

  double calc_b2() {
    double s = 0.0;
    for (int i = 1; i < root_no; ++i) {
      // we are only interested in tip depths
      s += tree[i].depth / std::pow(2, tree[i].depth);
    }
    return s;
  }

  double var_leaf_depth() {
    double average_depth = 0.0;
    int n = root_no - 1;
    for (int i = 1; i < root_no; ++i) {
      average_depth += tree[i].depth;
    }

    average_depth *= 1.0 / n;

    double var_depth = 0.0;
    for (int i = 1; i < root_no; ++i) {
      var_depth += (tree[i].depth - average_depth) *
        (tree[i].depth - average_depth);
    }
    var_depth *= 1.0 / n;
    return var_depth;
  }

  double calc_tot_int_path() {
    double sum_depth = 0.0;
    for (size_t i = static_cast<size_t>(root_no); i < tree.size(); ++i) {
      sum_depth += tree[i].depth;
    }
    return sum_depth;
  }

  double calc_avg_vert_depth() {
    double sum_depth = 0.0;
    for (size_t i = 1; i < tree.size(); ++i) {
      sum_depth += tree[i].depth;
    }
    auto answ =  sum_depth * 1.0 / (tree.size() - 1);
    return answ;
  }
};

}   // namespace width

