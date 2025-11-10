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

using ltable = std::vector< std::array<double, 4>>;

size_t find_daughters(const ltable& ltab_in,
                      double id,
                      double bt) {
  size_t num_daughters = 0;

  auto start = std::lower_bound(ltab_in.begin(), ltab_in.end(), bt,
                                [](const auto& a, double val) {
                                  return a[0] > val; });

  for (auto i = start; i != ltab_in.end(); ++i) {
    if ((*i)[1] == id && (*i)[0] <= bt) {
      num_daughters++;
      if (num_daughters > 1) break;
    }
  }
  return num_daughters;
}

size_t calc_cherries_ltable(const ltable& ltab_in) {
  size_t num_cherries = 0;
  for (const auto& i : ltab_in) {
    auto parent = i[1];

    if (parent == 0) continue;
    auto bt = i[0];

    size_t num_daughter_branches = find_daughters(ltab_in, parent, bt);
    size_t own_daughters         = find_daughters(ltab_in, i[2], bt);

    if (num_daughter_branches == 1 && own_daughters == 0) {
      num_cherries++;
    }
  }
  return num_cherries;
}
