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
#include <algorithm>
#include <vector>
#include "phylotree.h"   // NOLINT [build/include_subdir]

namespace max_depth {
class max_depth_tree {
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    size_t depth = 0;

    void set_depth() {
      if (!daughter1 && !daughter2) {
        depth = 1;
      } else {
        if (daughter1 && !daughter2) {
          depth = 1 + daughter1->depth;
        } else {
          auto d1 = daughter1->depth;
          auto d2 = daughter2->depth;
          depth = 1 + std::max(d1, d2);
        }
      }
    }
  };

  phylo_tree_t<node_t> tree;

 public:
  explicit max_depth_tree(const std::vector< int >& tree_edge)
    : tree(make_phylo_tree<node_t, false>(tree_edge)) {
  }

  int max_depth() {
    size_t md = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      (*i).set_depth();
      if ((*i).depth > md) md = (*i).depth;
    }
    return md;
  }

  double avg_depth() {
    double md = 0.0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      (*i).set_depth();
      md += (*i).depth;
    }
    md *= 1.0 / tree.size();
    return md;
  }
};

}  // end namespace max_depth

