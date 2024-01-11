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

std::vector< double > branching_times(
    const std::vector< std::array< size_t, 2>>& edge,
    const std::vector<double>& edge_length,
    size_t Nnode,
    size_t n) {
  std::vector<double> xx(Nnode, 0.0);

  for (size_t i = 0; i < edge_length.size(); ++i) {
    if (edge[i][1] > n) {
      auto target_index = edge[i][1] - n - 1;  // e2[i] - n, -2 because -1 of R,
                                              // and -1 of n (R->CPP conversion)
      auto source_index = edge[i][0] - n - 1;  // e1[i] - n
      xx[target_index] = xx[source_index] + edge_length[i];
    }
  }

  auto edge_index = edge[edge_length.size() - 1][0] - n - 1;

  double depth = xx[edge_index] + edge_length.back();

  for (auto& i : xx) {
    i = depth - i;
  }
  return xx;
}
