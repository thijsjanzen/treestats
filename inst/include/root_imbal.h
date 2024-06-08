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
#include <algorithm>

using ltable = std::vector< std::array<double, 4>>;

double calc_root_imbal(const ltable& ltab) {
  std::array<double, 2> lin_cnt = {0, 0};
  for (const auto& i : ltab) {
    if (i[2] < 0) {
      lin_cnt[0]++;
    } else {
      lin_cnt[1]++;
    }
  }
  double div = lin_cnt[0] + lin_cnt[1];
  double answ = lin_cnt[0] / div;
  if (answ < 0.5) answ = 1 - answ;

  return answ;
}


