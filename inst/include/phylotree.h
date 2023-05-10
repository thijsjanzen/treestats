// Copyright 2022 - 2023 Thijs Janzen
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
template <typename NODE>
auto make_phylo_tree(const std::vector<int>& tree_edge,
                     bool full_tree = false) {
  // this holds always:
  int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
  int tree_size = tree_edge.size() / 2 - root_no + 2;
  if (full_tree) tree_size = 1 + *std::max_element(tree_edge.begin(),
                                                   tree_edge.end());
  auto tree = phylo_tree_t<NODE>(tree_size);

  for (size_t i = 0; i < tree_edge.size(); i += 2) {
    int index    = static_cast<int>(tree_edge[i])     - root_no;
    int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

    if (full_tree) {
      index += root_no;
      d1_index += root_no;
    }

    if (full_tree || d1_index > 0) {
      !tree[index].daughter1 ?
       tree[index].daughter1 = &tree[d1_index] :
       tree[index].daughter2 = &tree[d1_index];
    }
  }
  return tree;
}
