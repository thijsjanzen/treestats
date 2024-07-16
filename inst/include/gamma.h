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
    total += (i + 1) * (brts_[i] - brts_[i - 1]);

    double_sum += total;
  }
  total += n * (h - brts_.back());

  double a = double_sum * 1.0 / (n - 2);

  double b = total * 0.5;
  double c = total * sqrt(1.f / (12 * n - 24));

  return (a - b) / c;
}

double calc_gamma2(const std::vector<int>& t_edge,
                   const std::vector<double>& edge_length) {
  int Nnode = static_cast<int>(edge_length.size()) / 2;
  int n = Nnode + 1;
  int n_1 = n + 1;

  std::vector<double> xx(Nnode, 0.f);

  for (size_t i = 0; i < t_edge.size(); i += 2) {
    auto j = i / 2;
    if (t_edge[i + 1] > n) {
      xx[t_edge[i + 1] - n_1] = xx[t_edge[i] - n_1] + edge_length[j];
    }
  }

  auto h = xx[t_edge[t_edge.size() - 2] - n_1] + edge_length.back();

  std::partial_sort(xx.begin(), xx.begin() + n - 1, xx.end());

  double total = 0.0;
  double double_sum = 0.0;

  for (int i = 1; i < n - 1; ++i) {
    total += (i + 1) * (xx[i] - xx[i - 1]);
    double_sum += total;
  }

  total += n * (h - xx[n - 2]);

  double c = total * sqrt(1.0 / (12 * n - 24));

  return (double_sum / (n - 2) - 0.5 * total) / c;
}
