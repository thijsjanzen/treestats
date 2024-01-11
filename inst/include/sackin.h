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
#include <algorithm>
#include <array>
#include <string>

#include "binom.h"        // NOLINT [build/include_subdir]
#include "phylotree.h"   // NOLINT [build/include_subdir]

using ltable = std::vector< std::array<double, 4>>;

double calc_sackin(const ltable& ltable_,
                   std::string normalization) {
  std::vector< int > s_values(ltable_.size(), 0);
  s_values[0] = 1;
  s_values[1] = 1;

  // ltable:
  // 0 = branching time // not used here
  // 1 = parent
  // 2 = id
  // 3 = extinct time // not used here
  for (size_t i = 2; i < ltable_.size(); ++i) {
    int parent_index = abs(static_cast<int>(ltable_[i][1])) - 1;
    s_values[parent_index]++;
    s_values[i] = s_values[parent_index];
  }
  // verified with R for correct values
  // compared with apTreeShape results, based on ltable
  // 24-09-2021
  double s = std::accumulate(s_values.begin(), s_values.end(), 0);

  if (normalization == "yule") {
    double sum_count = 0.0;
    size_t n = ltable_.size();
    for (size_t j = 2; j <= n; ++j) {
      sum_count += 1.0 / j;
    }
    return 1.0 * (s - 2.0 * n * sum_count) / n;
  }

  if (normalization == "pda") {
    size_t n = ltable_.size();
    double denom = powf(n, 1.5f);
    return 1.0 * s /denom;
  }

  return s;
}

namespace sackin {
class sackin_tree {
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    size_t num_extant_tips = 0;

    size_t get_acc_num_tips() {
      if (!daughter1 && !daughter2) {
        num_extant_tips = 2;
      } else {
        if (daughter1 && !daughter2) {
          num_extant_tips = 1 + daughter1->num_extant_tips;
        } else {
          num_extant_tips = daughter1->num_extant_tips +
            daughter2->num_extant_tips;
        }
      }
      return num_extant_tips;
    }
  };

  phylo_tree_t<node_t> tree;

 public:
  explicit sackin_tree(const std::vector< int >& tree_edge)
    : tree(make_phylo_tree<node_t, false>(tree_edge)) {
  }

  double calc_sackin() {
    double s = 0.0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      s += (*i).get_acc_num_tips();
    }
    return s;
  }

  double calc_tot_coph() {
    double tot_coph = 0.0;
    for (size_t i = tree.size() - 1; i >=  1; --i) {
      tree[i].get_acc_num_tips();
      if (tree[i].num_extant_tips > 0) {
        tot_coph += binom_coeff(tree[i].num_extant_tips, 2);
      }
    }
    return tot_coph;
  }

  size_t count_pitchforks() {
    size_t num_pitchforks = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      if ((*i).get_acc_num_tips() == 3) {
        num_pitchforks++;
      }
    }
    return num_pitchforks;
  }

  size_t count_cherries() {
    size_t num_cherries = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      if ((*i).get_acc_num_tips() == 2) {
        num_cherries++;
      }
    }
    return num_cherries;
  }

  double calc_blum() {
    double s = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      if ((*i).get_acc_num_tips() > 1) {
        s += log(1.0 * (*i).num_extant_tips - 1);
      }
    }
    return s;
  }
};

}  // end namespace sackin

namespace correction {
double correct_pda(size_t n,
                   double Is) {
  double denom = powf(n, 1.5f);
  return 1.0 * Is / denom;
}

double correct_yule(size_t n,
                    double Is) {
  double sum_count = 0.0;
  for (size_t j = 2; j <= n; ++j) {
    sum_count += 1.0 / j;
  }
  return 1.0 * (Is - 2.0 * n * sum_count) / n;
}

double correct_blum(size_t n,
                    double Is) {
  return Is * 1.0 / n;
}
}   // end namespace correction
