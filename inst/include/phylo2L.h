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

#include <Rcpp.h>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

#include <utility>

std::vector< double > branching_times_cpp(const Rcpp::List& phy);

size_t get_min_index(const std::vector< std::array<double, 6>>& localtab,
                     size_t col_index) {
  auto min_entry = std::min_element(localtab.cbegin(), localtab.cend(),
                                    [&](const auto& a, const auto& b) {
                                      return a[col_index] < b[col_index];
                                    });
  return std::distance(localtab.cbegin(), min_entry);
}

bool parent_in_nodesindex(const std::vector< size_t >& nodesindex,
                          size_t parent) {
  return std::binary_search(nodesindex.begin(), nodesindex.end(), parent);
}


void remove_from_L(std::vector< std::array<double, 6>>* L,
                   size_t j) {
  std::swap((*L)[j], (*L).back());
  (*L).pop_back();
}


std::vector< std::array<double, 6>> get_realL(
    const std::vector< size_t >& nodesindex,
    std::vector< std::array<double, 6>> L) {

  std::vector< std::array<double, 6>> realL;

  while (true) {
    size_t j = get_min_index(L, 2);
    size_t daughter = L[j][2];
    size_t parent   = L[j][1];
    if (parent_in_nodesindex(nodesindex, parent)) {
      for (auto& index : L) {
        if (index[1] == parent) {
          index[1] = daughter;
        }
      }

      bool match_found = false;
      for (auto& i : L) {
        if (i[2] == parent) {   // can only have one parent.
          i[5] = L[j][5];
          i[2] = daughter;
          remove_from_L(&L, j);
          match_found = true;
          break;
        }
      }
      if (!match_found) {
        realL.push_back(L[j]);
        remove_from_L(&L, j);
      }

    } else {
      realL.push_back(L[j]);
      remove_from_L(&L, j);
    }

    if (L.empty()) {
      break;
    }
  }

  std::sort(realL.begin(), realL.end(), [&](const std::array< double, 6>& v1,
                        const std::array< double, 6>& v2) {
      return v1[0] == v2[0] ? v1[2] < v2[2] : v1[0] > v2[0];
    });
  return realL;
}

size_t find_index(const std::vector< std::array<double, 6>>& pre_Ltable,
                  double ref) {
  size_t index = 0;
  for (; index != pre_Ltable.size(); ++index) {
    if (pre_Ltable[index][2] == ref) {
      break;
    }
  }
  return index;
}

std::vector< std::array< double, 4> > phylo_to_l_cpp(const Rcpp::List& phy) {
  std::vector< double > brts = branching_times_cpp(phy);

  auto min_brts = *std::min_element(brts.begin(), brts.end());
  if (min_brts < 0.0) {
    for (auto& i : brts) {
      i += fabs(min_brts);
    }
  }

  size_t num_species = static_cast<size_t>(phy["Nnode"]) + 1;

  Rcpp::StringVector tiplabel = phy["tip.label"];
  Rcpp::NumericMatrix edge    = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  size_t num_tips = tiplabel.size();

  std::vector< double > brt_preL(edge.nrow());
  double min_brt_preL = 1e10;

  for (int i = 0; i < edge.nrow(); ++i) {
    auto index = edge(i, 0) - num_tips - 1;   // -1 because 0 indexing
    brt_preL[i] = brts[index];
    if (brt_preL[i] < min_brt_preL) {
      min_brt_preL = brt_preL[i];
    }
  }

  if (min_brt_preL == 0.0) {
    double correction = 0.0;
    for (int i = 0; i < edge_length.size(); ++i) {
      if (brt_preL[i] == 0.0) {
        if (edge_length[i] > correction) {
          correction = edge_length[i];
        }
      }
    }

    for (auto& i : brt_preL) {
      i += correction;
    }
  }

  std::vector< std::array<double, 6>> pre_Ltable(brt_preL.size());

  for (size_t i = 0; i < brt_preL.size(); ++i) {
    pre_Ltable[i][0] = brt_preL[i];
    pre_Ltable[i][1] = edge(i, 0);
    pre_Ltable[i][2] = edge(i, 1);
    pre_Ltable[i][3] = edge_length[i];
    pre_Ltable[i][4] = brt_preL[i] - edge_length[i];
  }

  // all identical up to here (23-11-2023)
  std::vector<double> eeindicator(edge_length.size(), 0);

  std::vector< size_t > extant_species_index;
  for (size_t i = 0; i < pre_Ltable.size(); ++i) {
    if (pre_Ltable[i][4] <= 1e-6) {
      extant_species_index.push_back(pre_Ltable[i][2]);

      size_t index = find_index(pre_Ltable, pre_Ltable[i][2]);

      if (index < pre_Ltable.size()) {
        eeindicator[index] = -1;
      }
    }
  }

  std::sort(extant_species_index.begin(), extant_species_index.end());
  for (size_t i = 1; i <= num_species; ++i) {
    bool found = std::binary_search(extant_species_index.begin(),
                                    extant_species_index.end(),
                                    i);
    if (!found) {
      size_t index = find_index(pre_Ltable, i);

      if (index < pre_Ltable.size()) {
        eeindicator[index] = pre_Ltable[index][4];
      }
    }
  }

  for (const auto& i : extant_species_index) {
    size_t index = find_index(pre_Ltable, i);

    if (index < pre_Ltable.size()) {
      eeindicator[index] = -1;
    }
  }

  // verified correct so far
  for (size_t i = 0; i < eeindicator.size(); ++i) {
    pre_Ltable[i][5] = eeindicator[i];
  }

  // verified identical

  std::sort(pre_Ltable.begin(), pre_Ltable.end(),
            [&](const std::array< double, 6>& v1,
                const std::array< double, 6>& v2) {
    return(v1[0] > v2[0]);   // sort decreasing
  });

  std::vector< size_t > nodesindex(edge.nrow());
  for (int i = 0; i < edge.nrow(); ++i) {
    nodesindex[i] = static_cast<size_t>(edge(i, 0));
  }

  std::sort(nodesindex.begin(), nodesindex.end());
  nodesindex.erase(std::unique(nodesindex.begin(), nodesindex.end()),
                   nodesindex.end());

  std::vector< std::array<double, 6>> realL = get_realL(nodesindex,
                                                        pre_Ltable);

  std::vector< std::array< double, 4> > L(realL.size());

  for (size_t i = 0; i < realL.size(); ++i) {
    L[i][0] = realL[i][0];
    L[i][1] = realL[i][1];
    L[i][2] = i + 1;
    L[i][3] = realL[i][5];

    size_t index = 0;
    for (; index < realL.size(); ++index) {
      if (realL[index][2] == realL[i][1]) {
        break;
      }
    }
    if (index != realL.size()) {
      L[i][1] = index + 1;
    }
  }

  L[0][1] = 0;
  L[0][2] = -1;
  L[1][1] = -1;

  for (size_t i = 1; i < L.size(); ++i) {
    if (L[i - 1][2] < 0) {
      auto ref = std::abs(L[i - 1][2]);
      for (auto& j : L) {
        if (j[1] == ref) {
          j[1] = L[i - 1][2];
          j[2] *= -1;
        }
      }
    }
  }

  return L;
}
