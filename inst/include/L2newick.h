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

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <utility>

size_t which_max_index(const std::vector< std::array< double, 4>>& ltable) {
  auto max_val = std::max_element(ltable.begin(), ltable.end(),
                                  [&](const auto& a, const auto& b) {
                                      return a[0] < b[0];});
  return std::distance(ltable.begin(), max_val);
}

int index_of_parent(const std::vector< std::array< double, 4>>& ltable,
                    int parent) {
  int index = 0;
  bool found = false;
  for (; index < static_cast<int>(ltable.size()); ++index) {
    if (std::abs(ltable[index][2] - parent) < 0.0000001) {
      found = true;
      break;
    }
  }
  if (!found) {
    index = -1;
  }
  return index;
}

std::string d_to_s(double d) {
  std::stringstream out;
  out << std::fixed << std::setprecision(15) << d;
  return out.str();
}

void remove_from_dataset(std::vector< std::array< double, 4>>* ltable,
                         std::vector< std::string>* linlist,
                         size_t index) {
  std::swap((*ltable)[index], (*ltable).back());
  (*ltable).pop_back();
  std::swap((*linlist)[index], (*linlist).back());
  (*linlist).pop_back();
}

std::string ltable_to_newick(const std::vector< std::array< double, 4>>& ltable,
                       bool drop_extinct) {
  auto L = ltable;
  //  first sort ltable
  //  L = L[order(abs(L[, 3])), 1:4]
  std::sort(L.begin(), L.end(), [&](const auto& a, const auto& b) {
    return std::abs(static_cast<int>(a[2])) < std::abs(static_cast<int>(b[2]));
  });

  double age = L[0][0];   // age = L[1, 1]
  std::vector< std::array< double, 4>> new_L;

  for (auto& i : L) {
    i[0] = age - i[0];   // L[, 1] = age - L[, 1]
    if (i[0] < 0.0) i[0] = 0.0;

    // i[3] == -1
    bool is_extant = (std::abs(i[3] + 1) < 0.0000001);

    if (!is_extant) {   // notmin1 = which(L[, 4] != -1)
      // L[notmin1, 4] = age - L[notmin1, 4]
      i[3] = age - i[3];
    } else {
      i[3] = age;
    }

    if (drop_extinct == true) {
      if (is_extant) {
        new_L.push_back(i);
      }
    }
  }

  std::vector< std::array< double, 4>> L_original = L;

  if (drop_extinct == true) {
      // keep a copy of the original ltable for later look up purpose
      //   L_original = L;
    L = new_L;
    L_original[0][0] = -1.0;
  }

  // check whether L[1, ] is t1
  // if yes, change L[1, 1] to -1
  if (std::abs(L[0][1]) < 0.0000001) {
    L[0][0] = -1.0;
  }

  std::vector< std::string > linlist_4(L.size());
  size_t index = 0;
  for (const auto& i : L) {
    std::string add = "t" + std::to_string(abs(static_cast<int>(i[2])));
    linlist_4[index] = add;
    index++;
  }

  if (linlist_4.size() != L.size()) {
    throw std::invalid_argument("linlist_4.size() != L.size()");
  }

  // verified correct for 4/5 tip phylogeny up until here.
  while (true) {
    auto j = which_max_index(L);
    int parent    = static_cast<int>(L[j][1]);
    int parentj   = index_of_parent(L, parent);
    if (parentj != -1) {   // -1 means not found
      double bl = std::abs(L[parentj][3] - L[j][0]);
      std::string spec1 = linlist_4[parentj] + ":" + d_to_s(bl);
      double bl2 = std::abs(L[j][3] - L[j][0]);
      std::string spec2 = linlist_4[j] + ":" + d_to_s(bl2);
      linlist_4[parentj] = "(" + spec1 + "," + spec2 + ")";
      L[parentj][3] = L[j][0];
      remove_from_dataset(&L, &linlist_4, j);
    } else {
      parentj = index_of_parent(L_original, parent);
      if (parentj == -1) {
        throw std::invalid_argument("Look up failed "+ std::to_string(j) +
                                     " " + std::to_string(parent) + "\n");
      }
      for (int i = 0; i <= 2; ++i) {
        L[j][i] = L_original[parentj][i];
      }
    }

    if (linlist_4.size() == 1) {
      break;
    }
  }

  return linlist_4[0] + ":" + d_to_s(L[0][3]) + ";";
}
