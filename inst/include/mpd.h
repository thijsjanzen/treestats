// Copyright 2022 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#pragma once
#include <vector>

namespace mpd_tree {

// here, we don't use phylo_tree.h, because we also need to
// track edge lengths. Code is very similar though.
struct mpd_node {
  mpd_node* daughterL = nullptr;
  mpd_node* daughterR = nullptr;
  size_t L;
  size_t R;
  double bl_R;
  double bl_L;

  mpd_node() {
    daughterL = nullptr;
    daughterR = nullptr;
    L = 1;
    R = 1;
    bl_R = -1.0;
    bl_L = -1.0;
  }

  size_t update_num_tips() {
    if (daughterL && !daughterR) {
      L = daughterL->update_num_tips();
    }
    if (daughterR && !daughterL) {
      R = daughterR->update_num_tips();
    }
    if (daughterL && daughterR) {
      L =  daughterL->update_num_tips();
      R =  daughterR->update_num_tips();
    }

    return L + R;
  }
};

class phylo_tree {
 public:
  explicit phylo_tree(const std::vector< int >& tree_edge,
                      const std::vector<double>& edge_length) {
    int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
    tree_size = root_no - 1;

    tree.resize(tree_edge.size() / 2 - root_no + 2);

    for (size_t i = 0; i < tree_edge.size(); i += 2) {
      int index    = static_cast<int>(tree_edge[i]) - root_no;
      int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;
      int el_index = i / 2;

      if (d1_index >= 0) {
        // we are dealing with an internal node
        tree[index].bl_L < 0 ?  // was the left daughter already set before?
          tree[index].daughterL = &tree[d1_index] :
          tree[index].daughterR = &tree[d1_index];
      }

      tree[index].bl_L < 0 ?
        tree[index].bl_L = edge_length[el_index] :
        tree[index].bl_R = edge_length[el_index];
    }

    tree[0].update_num_tips();
  }

  double calculate_mpd() {
    int N = tree_size;
    double mpd = 0.0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      auto L_bl = i.bl_L;
      auto R_bl = i.bl_R;

      mpd += L_bl * (l * (N - l));
      mpd += R_bl * (r * (N - r));
    }

    mpd *= 2.0 / (N * (N - 1));
    return(mpd);
  }

 private:
  std::vector< mpd_node > tree;
  int tree_size = 0;
};


}  // namespace mpd_tree
