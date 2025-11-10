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
#include <numeric>
#include <algorithm>  // std::sort

const double prefactor = 2 * sqrt(3);

double calc_gamma(std::vector<double> brts_) {
  double n = brts_.size() + 1;

  auto h = *std::max_element(brts_.begin(), brts_.end());

  for (auto& i : brts_) {
    i =  h - i;
  }

  std::sort(brts_.begin(), brts_.end());

  double total = 0.0;
  double double_sum = 0.0;
  auto max_i = n - 1;
  for (size_t i = 1; i < max_i; ++i) {
    if (i - 1 < 0) throw "gamma: out of range, i-1 < 0";
    if (i >= brts_.size()) throw "gamma: out of range, i >= brts_.size()";

    total += (i + 1) * (brts_[i] - brts_[i - 1]);

    double_sum += total;
  }
  total += n * (h - brts_.back());

  double a = double_sum * 1.0 / (n - 2);

  double b = total * 0.5;
  double c = total * sqrtf(1.f / (12 * n - 24));

  return (a - b) / c;
}
