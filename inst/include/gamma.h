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
#include <numeric>
#include <algorithm>  // std::sort

const double prefactor = 2 * sqrt(3);

double calc_gamma(std::vector<double> brts_) {
  double n = brts_.size() + 1;

  auto h = brts_[0];   // *std::max_element(brts_.begin(), brts_.end());

  for (auto& i : brts_) {
    i =  h - i;
  }

  std::sort(brts_.begin(), brts_.end());

  double total = 0.0;
  double double_sum = 0.0;
  double temp = n * (h - brts_.back());
  std::adjacent_difference(brts_.begin(), brts_.end(), brts_.begin());

  size_t j = 1;
  for (const auto& i : brts_) {
    total += j * i;
    double_sum += total;
    j++;
  }

  total += temp;

  double mult_total = 1.0 / total;
  return prefactor *
          sqrt(n - 2) *
          (double_sum * 1.0 / (n - 2) - total * 0.5) *
          mult_total;
}
