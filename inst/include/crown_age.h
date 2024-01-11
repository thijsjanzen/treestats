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

#include <functional>  // std::greater
#include <algorithm>   // std::sort
#include <vector>

#include "phylodiv.h"  // NOLINT [build/include_subdir]


double get_total_bl(const std::vector< std::array<size_t, 2>>& edge,
                    const std::vector<double>& el,
                    size_t init_label) {
  size_t tip_label = init_label;
  size_t root_label = edge[0][0];
  size_t tip_index = 0;
  for (tip_index = 0; tip_index < edge.size(); ++tip_index) {
    if (edge[tip_index][1] == tip_label) {
      break;
    }
  }
  double bl = el[tip_index];
  while (edge[tip_index][0] != root_label) {
    tip_label = edge[tip_index][0];
    for (tip_index = 0; tip_index < edge.size(); ++tip_index) {
      if (edge[tip_index][1] == tip_label) {
        break;
      }
    }
    bl += el[tip_index];
  }
  return bl;
}

void update_dist_to_root(std::vector<double>* dist_to_root,
                         size_t* focal_index,
                         const std::vector< std::array<size_t, 2>>& edge,
                         const std::vector<double>& el) {
  double bl = get_total_bl(edge, el, *focal_index);
  (*focal_index)++;
  (*dist_to_root).push_back(bl);
  std::sort((*dist_to_root).begin(), (*dist_to_root).end(),
            std::greater<double>());
  return;
}

double calc_crown_age(std::vector< std::array<size_t, 2>> edge,
                      std::vector<double> el) {
  sort_edge_and_edgelength(&edge, &el);

  size_t focal_index = 1;
  size_t root_label = edge[0][0];

  std::vector<double> dist_to_root;
  update_dist_to_root(&dist_to_root, &focal_index, edge, el);
  update_dist_to_root(&dist_to_root, &focal_index, edge, el);

  while (dist_to_root[1] != dist_to_root[0] && focal_index < root_label) {
    update_dist_to_root(&dist_to_root, &focal_index, edge, el);
  }
  return dist_to_root[0];
}
