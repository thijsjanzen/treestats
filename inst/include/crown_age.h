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

double calc_crown_age(std::vector< std::array<size_t, 2>> edge,
                      std::vector<double> el) {
  // we calculate the distance to the root first for all nodes
  // then for all tips
  size_t max_label_id = 1 + el.size();

  //  + 1 because of R indexing
  std::vector<double> dist_to_root(max_label_id + 1, 0.0);
  double crown_age = -1;
  for (size_t i = 0; i < edge.size(); ++i) {
    size_t focal_index  = edge[i][1];
    size_t parent_index = edge[i][0];
    double bl = el[i];

    auto d = bl + dist_to_root[parent_index];

    dist_to_root[focal_index] = d;

    if (d > crown_age) crown_age = d;
  }

  return crown_age;
}
