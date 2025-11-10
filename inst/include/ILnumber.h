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
//
#pragma once
#include <vector>

size_t calc_IL(const std::vector< int >& tree_edge) {
  int root_no = 2 + tree_edge.size() * 0.25;
  std::vector< int > nodes(tree_edge.size() / 2, 0);

  for (size_t i = 0; i < tree_edge.size(); i += 2) {
    if (tree_edge[i + 1] < root_no) {
      int index    = static_cast<int>(tree_edge[i]) - root_no;
      nodes[index]++;
    }
  }

  size_t num_IL = 0;
  for (const auto& i : nodes) {
    if (i == 1) num_IL++;
  }
  return(num_IL);
}
