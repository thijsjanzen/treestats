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
#include <algorithm>
#include <vector>

#include "phylotree.h"  // NOLINT [build/include_subdir]

namespace b1_tree {
class b1_tree {
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    int max_dist_to_tips = 1;

    void update_max_dist_to_tips() {
      if (daughter1 && !daughter2) {
        max_dist_to_tips = 1 + daughter1->max_dist_to_tips;
      } else if (daughter1 && daughter2) {
        auto d1 = 1 + daughter1->max_dist_to_tips;
        auto d2 = 1 + daughter2->max_dist_to_tips;
        max_dist_to_tips = std::max(d1, d2);
      }
      return;
    }
  };

  phylo_tree_t<node_t> tree;

 public:
  explicit b1_tree(const std::vector< int >& tree_edge)
    : tree(make_phylo_tree<node_t, false>(tree_edge)) {
  }

  double calc_b1() {
    double b1 = 0.0;
    for (size_t i = tree.size() - 1; i >= 1; --i) {
      tree[i].update_max_dist_to_tips();
      b1 += 1.0 / tree[i].max_dist_to_tips;
    }
    return b1;
  }
};

}   // namespace b1_tree
