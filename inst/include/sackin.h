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
#include <algorithm>
#include <array>
#include <string>

#include "binom.h" // NOLINT [build/include_subdir]

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

struct node {
  node* daughter1 = nullptr;
  node* daughter2 = nullptr;
  size_t num_extant_tips;

  node() {
    num_extant_tips = 0;
  }

  size_t get_acc_num_tips() {
    if (!daughter1 && !daughter2) {
      num_extant_tips = 2;
    } else {
      if (daughter1 && !daughter2) {
        num_extant_tips = 1 + daughter1->get_acc_num_tips();
      } else {
        num_extant_tips = daughter1->get_acc_num_tips() +
                          daughter2->get_acc_num_tips();
      }
    }
    return num_extant_tips;
  }
};

class phylo_tree {
 public:
  explicit phylo_tree(const std::vector< int >& tree_edge) {
    int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
    tree.resize(tree_edge.size() / 2);

    for (size_t i = 0; i < tree_edge.size(); i += 2) {
      int index    = static_cast<int>(tree_edge[i]) - root_no;
      int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

      if (d1_index > 0) {
        !tree[index].daughter1 ?
         tree[index].daughter1 = &tree[d1_index] :
         tree[index].daughter2 = &tree[d1_index];
      }
    }
  }

  int calc_sackin() {
    tree[0].get_acc_num_tips();
    int s = 0;
    for (const auto& i : tree) {
      s += i.num_extant_tips;
    }
    return s;
  }

  double calc_tot_coph() {
    tree[0].get_acc_num_tips();
    double tot_coph = 0.0;
    for (size_t i = 1; i < tree.size(); ++i) {
      if (tree[i].num_extant_tips > 0) {
        tot_coph += binom_coeff(tree[i].num_extant_tips, 2);
      }
    }
    return tot_coph;
  }

  size_t count_pitchforks() {
    tree[0].get_acc_num_tips();
    size_t num_pitchforks = 0;
    for (const auto& i : tree) {
      if (i.num_extant_tips == 3) {
        num_pitchforks++;
      }
    }
    return num_pitchforks;
  }

  size_t count_cherries() {
    tree[0].get_acc_num_tips();
    size_t num_cherries = 0;
    for (const auto& i : tree) {
      if (i.num_extant_tips == 2) {
        num_cherries++;
      }
    }
    return num_cherries;
  }

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

  double calc_blum(bool normalize, size_t n) {
    tree[0].get_acc_num_tips();
    double s = 0;
    for (const auto& i : tree) {
      if (i.num_extant_tips > 1) {
        s += log(1.0 * i.num_extant_tips - 1);
      }
    }
    if (normalize) {
      s *= 1.0 / n;
    }
    return s;
  }

 private:
  std::vector< node > tree;
};
