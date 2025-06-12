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
#include <algorithm>

std::vector<double> create_normalized_brts(const std::vector<double>& v) {
  std::vector<double> output = v;
  std::sort(output.begin(), output.end());
  if (output.front() == 0.f) {
    // t goes from 0 to T
    double maxT = output.back();
    output.push_back(maxT);
    for (auto& i : output) {
      i = i / maxT;
    }
  } else {
    output.push_back(0.f);
    // we assume here that the brts go from -T to 0 ...
    auto min_brts = output.front();
    for (auto& i : output) {
      i = 1.f - i / min_brts;
    }
  }
  return output;
}

std::vector< double > create_normalized_lins(size_t num_lin) {
  std::vector< double > output(num_lin - 1);
  std::iota(output.begin(), output.end(), 2.f);
  auto max_num = output.back();
  output.push_back(max_num);
  for (auto& i : output) {
    i *= 1.f / max_num;
  }
  return(output);
}

int get_index(const std::vector<double>& brts,
              double tim) {
  auto index = std::lower_bound(brts.begin(), brts.end(), tim);
  if (index != brts.begin()) index--;

  return std::distance(brts.begin(), index);
}

double calc_nltt_from_data(const std::vector<double>& b1,
                           const std::vector<double>& b2,
                           const std::vector<double>& n1,
                           const std::vector<double>& n2,
                           const std::vector<double>& all_b) {
  double nltt = 0.f;
  for (size_t k = 1; k < all_b.size(); ++k) {
    double tim = all_b[k];
    auto index1 = get_index(b1, tim);
    auto index2 = get_index(b2, tim);

    auto num_lin1 = n1[index1];
    auto num_lin2 = n2[index2];
    auto diff_lin = num_lin1 - num_lin2;

    if (diff_lin < 0) diff_lin *= -1;
    double dt = all_b[k] - all_b[k-1];
    nltt += dt * diff_lin;
  }
  return nltt;
}

// please note that the branching times have to be from -T to 0
double calc_nltt(const std::vector<double>& v1,
                 const std::vector<double>& v2) {
  std::vector< double > b_times_1 = create_normalized_brts(v1);
  std::vector< double > b_times_2 = create_normalized_brts(v2);

  std::vector< double > lin_1 = create_normalized_lins(v1.size());
  std::vector< double > lin_2 = create_normalized_lins(v2.size());

  std::vector< double > all_branching_times(b_times_1.size() +
                                            b_times_2.size());
  std::merge(b_times_1.begin(), b_times_1.end(),
             b_times_2.begin(), b_times_2.end(),
             all_branching_times.begin());
  // notice that this might introduce duplicate branching times,
  // but this does not
  // affect nltt, as dt = 0.
  return calc_nltt_from_data(b_times_1, b_times_2,
                             lin_1, lin_2,
                             all_branching_times);
}
