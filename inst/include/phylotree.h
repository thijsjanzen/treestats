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

template <typename NODE>
using phylo_tree_t = std::vector<NODE>;

// standalone 'virtual constructor'
template <typename NODE,
          bool FULL_TREE>
auto make_phylo_tree(const std::vector<int>& tree_edge) {
  int root_no = tree_edge[0];
  for (size_t i = 2; i < tree_edge.size(); i+=2) {
    if (tree_edge[i] < root_no) root_no = tree_edge[i];
  }

  int tree_size = tree_edge.size() / 2 - root_no + 2;

  if constexpr (FULL_TREE) tree_size = 2 + 0.5 * tree_edge.size();

  auto tree = phylo_tree_t<NODE>(tree_size);

  for (size_t i = 0; i < tree_edge.size(); i += 2) {
    int index, d1_index;

    if constexpr (FULL_TREE) {
      index = static_cast<int>(tree_edge[i]);
      d1_index = static_cast<int>(tree_edge[i + 1]);
    } else {
      index = static_cast<int>(tree_edge[i])     - root_no;
      d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;
    }

    if constexpr (FULL_TREE) {
      !tree[index].daughter1 ?
      tree[index].daughter1 = &tree[d1_index] :
      tree[index].daughter2 = &tree[d1_index];
    } else {
      if (d1_index > 0) {
        !tree[index].daughter1 ?
        tree[index].daughter1 = &tree[d1_index] :
        tree[index].daughter2 = &tree[d1_index];
      }
    }
  }
  return tree;
}
