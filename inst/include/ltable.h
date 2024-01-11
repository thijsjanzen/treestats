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

#include "binom.h"   // NOLINT [build/include_subdir]

using ltable = std::vector< std::array<double, 4>>;

namespace ltab {

class stat {
 public:
  explicit stat(const ltable& ltab_in) : ltable_(ltab_in) {}

  size_t find_parent(const ltable& ltable_,
                     int focal_id,
                     int start_index) {
    for (int i = start_index; i >= 0; i--) {
      if (static_cast<int>(ltable_[i][2]) == focal_id) {
        return i;
      }
    }

    if (start_index != static_cast<int>(ltable_.size())) {
      return find_parent(ltable_, focal_id, ltable_.size());
    } else {
      throw "out of bounds";
    }
  }

  double calc_tot_coph() {
    std::vector< int > tips_tracker(ltable_.size(), 1);
    std::vector< int > node_tips;
    for (int focal_index = ltable_.size() - 1;
         focal_index > 1;
         --focal_index) {
      // parent index is in R form.
      auto parent = std::abs(ltable_[focal_index][1]) - 1;
      auto num_tips = tips_tracker[focal_index] +
        tips_tracker[parent];
      tips_tracker[parent] = num_tips;
      node_tips.push_back(num_tips);
    }

    double tot_coph = 0.0;
    for (size_t i = 0; i < node_tips.size(); ++i) {
      if (node_tips[i] > 0) {
        tot_coph += binom_coeff(node_tips[i], 2);
      }
    }
    return tot_coph;
  }

  double calc_blum(bool normalize) {
    std::vector< int > s_values(ltable_.size(), 1);
    for (size_t i = ltable_.size() - 1; i > 0; i--) {
      int parent_index = abs(static_cast<int>(ltable_[i][1])) - 1;
      s_values[parent_index] += s_values[i];
      s_values[i] = s_values[parent_index];
    }

    double s = 0.0;
    for (size_t i = 1; i < s_values.size(); ++i) {
      if (s_values[i] != 0.0) {
        s += log(1.0 * s_values[i] - 1.0);
      }
    }

    if (normalize) {
      s *= 1.0 / ltable_.size();
    }

    return s;
  }


  auto collect_depths() {
    std::vector< int > s_values(ltable_.size(), 0);
    s_values[0] = 1;
    s_values[1] = 1;

    for (size_t i = 2; i < ltable_.size(); ++i) {
      int parent_index = abs(static_cast<int>(ltable_[i][1])) - 1;
      s_values[parent_index]++;
      s_values[i] = s_values[parent_index];
    }
    return s_values;
  }

  size_t calc_max_depth() {
    auto s_values = collect_depths();
    return(*std::max_element(s_values.begin(), s_values.end()));
  }

  double calc_b2() {
    auto depths = collect_depths();
    double s = 0.0;
    for (const auto& i : depths) {
      s += i / std::pow(2, i);
    }
    return s;
  }

  double calc_b1() {
    std::vector< int > depth_tracker(ltable_.size(), 1);
    double b1 = 0.0;
    std::vector<int> depths;
    for (int focal_index = ltable_.size() - 1;
         focal_index > 1;
         --focal_index) {
      // parent index is in R form.
      auto parent = std::abs(ltable_[focal_index][1]) - 1;
      auto max_dist = std::max(depth_tracker[focal_index],
                               depth_tracker[parent]);
      depth_tracker[parent] = 1 + max_dist;
      b1 += 1.0 / max_dist;
    }

    return b1;
  }

  double calc_var_leaf_depth() {
    std::vector< int > depths = collect_depths();

    double inv_tree_size = 1.0 / depths.size();
    double average_depth =
      std::accumulate(depths.begin(), depths.end(), 0.0) * inv_tree_size;

    double var_depth = 0.0;
    for (const auto& i : depths) {
      var_depth += (i - average_depth) * (i - average_depth);
    }
    var_depth *= 1.0 / (depths.size() - 1);
    return var_depth;
  }

  std::vector< int > collect_widths() {
    std::vector< int > current_depths(ltable_.size() + 1, 0);
    for (size_t i = 1; i < ltable_.size(); ++i) {
      int parent  = std::abs(ltable_[i][1]);
      int self_id = std::abs(ltable_[i][2]);
      current_depths.push_back(current_depths[parent]);
      current_depths[parent]++;
      current_depths[self_id] = current_depths[parent];
    }

    // make histogram
    std::vector<int> counts(ltable_.size(), 0);
    for (const auto& i : current_depths) {
      counts[i]++;
    }
    return counts;
  }

  size_t calc_max_width() {
    auto widths = collect_widths();
    return(*std::max_element(widths.begin(), widths.end()));
  }

  size_t max_del_width() {
    auto widths = collect_widths();
    std::vector<int> dW(widths.size() - 1);
    for (size_t i = 1; i < widths.size(); ++i) {
      dW[i - 1] = widths[i] - widths[i - 1];
    }
    return(*std::max_element(dW.begin(), dW.end()));
  }

 private:
  const ltable ltable_;
};
}   // namespace ltab
