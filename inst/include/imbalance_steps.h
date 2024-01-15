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

#include "L2newick.h"  // NOLINT [build/include_subdir]

using ltable = std::vector< std::array<double, 4>>;

namespace imbal_steps {

int get_attractor(const ltable& ltab) {
  int attractor = 2;
  int cnt_clade_1 = 0;
  int cnt_clade_2 = 0;
  for (const auto& i : ltab) {
    if (i[1] == -1) {
      cnt_clade_1++;
    } else if (i[1] == 2) {
      cnt_clade_2++;
    }
  }

  if (cnt_clade_1 > cnt_clade_2) attractor = -1;

  return attractor;
}

bool all_identical(const std::array< int, 5>& a) {
  for (size_t i = 1; i < 5; ++i) {
    if (a[i] != a[i - 1]) return false;
  }
  return true;
}

std::vector<int> get_daughters(const ltable& ltab,
                               int parent,
                               int focal_index) {
  std::vector<int> out;
  //  daughters <- which(new_ltab[, 2] == parent &
  //.                    new_ltab[, 1] > new_ltab[focal_index, 1])
  for (size_t i = 0; i < ltab.size(); ++i) {
    if (ltab[i][1] == parent) {
      if (ltab[i][0] > ltab[focal_index][0]) {
        out.push_back(i);
      }
    }
  }
  return out;
}

std::vector<int> find_daughters(const ltable& ltab,
                                int current_label,
                                int focal_index) {
  std::vector<int> out;
  //  aughters <- which(ltab[, 2] == current_label &
  //                    ltab[, 1] <= ltab[i, 1])
  for (size_t i = 0; i < ltab.size(); ++i) {
    if (ltab[i][1] == current_label) {
      if (ltab[i][0] <= ltab[focal_index][0]) {
        out.push_back(i);
      }
    }
  }
  return out;
}

std::vector<int> find_others(const ltable& ltab,
                             int current_label,
                             int focal_index) {
  std::vector<int> out;
  //  daughters <- which(ltab[, 2] == current_label &
  //                    ltab[, 1] <= ltab[i, 1])
  for (size_t i = 0; i < ltab.size(); ++i) {
    if (ltab[i][1] == current_label) {
      if (ltab[i][0] < ltab[focal_index][0]) {
        out.push_back(i);
      }
    }
  }
  return out;
}

void renumber_ltable(ltable* ltab) {
  auto temp_new_ltab = *ltab;

  for (size_t i = 0; i < temp_new_ltab.size(); ++i) {
    auto current_label = (*ltab)[i][2];
    if (std::abs(current_label) != (i + 1)) {
      int new_label = i + 1;   // +1 to adhere to R counting
      if (current_label < 0) new_label *= -1;
      temp_new_ltab[i][2] = new_label;
      auto daughters = find_daughters((*ltab), current_label, i);
      if (!daughters.empty()) {
        for (const auto& j : daughters) {
          temp_new_ltab[j][1] = new_label;
        }
      }

      auto other_instances = find_others((*ltab), i + 1, i);
      if (!other_instances.empty()) {
        for (const auto& j : other_instances) {
          temp_new_ltab[j][1] = current_label;
        }
      }
    }
  }

  *ltab = temp_new_ltab;
  return;
}

ltable swap_deepest(const ltable& ltab, int* main_attractor, bool* stop) {
  std::vector<int> depths(ltab.size(), 0);
  depths[0] = depths[1] = 1;
  for (size_t i = 2; i < ltab.size(); ++i) {
    int parent_index = -1 + std::abs(ltab[i][1]);
    depths[parent_index]++;
    depths[i] = depths[parent_index];
  }

  int max_depth_lineage = 1 +
    std::distance(depths.begin(),
                  std::max_element(depths.begin(), depths.end()));

  *main_attractor = 0;
  int focal_index = 0;
  for (focal_index = 0; focal_index < ltab.size(); ++focal_index) {
    if (std::abs(ltab[focal_index][2]) == max_depth_lineage) {
      *main_attractor = ltab[focal_index][2];
      break;
    }
  }

  auto new_ltab = ltab;
  *stop = false;

  if (std::abs(*main_attractor) > 2) {
    auto parent = ltab[focal_index][1];
    auto index_parent = index_of_parent(ltab, parent);
    new_ltab[focal_index][2] = parent;
    new_ltab[focal_index][1] = *main_attractor;
    new_ltab[index_parent][2] = *main_attractor;

    auto daughters = get_daughters(new_ltab, parent, focal_index);
    if (!daughters.empty()) {
      for (const auto& i : daughters) {
        new_ltab[i][1] = *main_attractor;
      }
    }

    renumber_ltable(&new_ltab);

  } else {
    *stop = true;
  }

  return new_ltab;
}



void rebase_ltable(ltable* ltab) {
  if (ltab->size() == 2) return;

  std::array< int, 5> prev_main_attractor;

  int current_main_attractor = -1;
  bool stop = false;
  size_t cnt = 0;
  while (!stop) {
    *ltab = swap_deepest(*ltab, &current_main_attractor, &stop);
    prev_main_attractor[cnt % prev_main_attractor.size()] =
                                                    current_main_attractor;

    cnt++;
    if (cnt > 3 && all_identical(prev_main_attractor)) {
      throw "Stuck in endless loop, possibly due to polytomies";
    }
  }

  renumber_ltable(ltab);
  return;
}

int number_of_steps(ltable ltab, bool normalization) {
  rebase_ltable(&ltab);

  auto attractor = get_attractor(ltab);

  int cnt_steps = 0;
  for (size_t i = 2; i < ltab.size(); ++i) {
    if (ltab[i][1] != attractor) cnt_steps++;
  }

  if (normalization) {
    size_t tree_size = ltab.size();
    int max_expected = tree_size - std::ceil(std::log2(tree_size)) - 1;
    cnt_steps /= max_expected;
  }

  return cnt_steps;
}
}   // namespace imbal_steps
