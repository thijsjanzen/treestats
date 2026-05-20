// Copyright 2022 - 2025 Thijs Janzen
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
    // int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
    int root_no = tree_edge[0];
    for (size_t i = 2; i < tree_edge.size(); i+=2) {
      if (tree_edge[i] < root_no) root_no = tree_edge[i];
    }

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
    // to prevent integer overflow when N * (N - 1) is too large
    double Nd = static_cast<double>(N);
    mpd *= 2.0 / (Nd * (Nd - 1));
    return(mpd);
  }

 private:
  std::vector< mpd_node > tree;
  int tree_size = 0;
};


// the path node is highly similar to the mpd node,
// just stores less

struct path_node {
  path_node* daughterL = nullptr;
  path_node* daughterR = nullptr;
  double bl_R;
  double bl_L;
  int id_R;
  int id_L;
  double add_one;

  // for inv_path length
  double sum_inv_b_length;

  path_node() {
    daughterL = nullptr;
    daughterR = nullptr;
    bl_R = -1.0;
    bl_L = -1.0;
    sum_inv_b_length = 0.0;
    id_R = -1;
    id_L = -1;
    add_one = 1.0;
  }

  void update_path(double branch, double prev_sum) {
    sum_inv_b_length += branch > 0.0 ? prev_sum + 1.0 / (add_one + branch) :
                                       prev_sum;

    if (daughterL) {
      daughterL->update_path(bl_L, sum_inv_b_length);
    }
    if (daughterR) {
      daughterR->update_path(bl_R, sum_inv_b_length);
    }
  }

  std::vector<std::vector<double>> get_inv_blengths() const {
    std::vector<std::vector<double>> v;
    if (!daughterL) v.push_back({sum_inv_b_length + 1.0 / (add_one + bl_L),
                                 static_cast<double>(id_L)});
    if (!daughterR) v.push_back({sum_inv_b_length + 1.0 / (add_one + bl_R),
                                 static_cast<double>(id_R)});

    return v;
  }
};

class phylo_path_tree {
 public:
  explicit phylo_path_tree(const std::vector< int >& tree_edge,
                           const std::vector<double>& edge_length,
                           bool add_one) {
    // int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
    int root_no = tree_edge[0];
    for (size_t i = 2; i < tree_edge.size(); i+=2) {
      if (tree_edge[i] < root_no) root_no = tree_edge[i];
    }

    tree_size = root_no - 1;

    tree.resize(tree_edge.size() / 2 - root_no + 2);
    for (auto& i : tree) {
        i.add_one = static_cast<double>(add_one);
    }

    for (size_t i = 0; i < tree_edge.size(); i += 2) {
      int index    = static_cast<int>(tree_edge[i]) - root_no;
      int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;
      int el_index = i / 2;

      if (d1_index >= 0) {
        // we are dealing with an internal node
        tree[index].bl_L < 0 ?  // was the left daughter already set before?
          tree[index].daughterL = &tree[d1_index] :
          tree[index].daughterR = &tree[d1_index];
      } else {
        // external node, store the id of the tip
        tree[index].bl_L < 0 ?
          tree[index].id_L = tree_edge[i + 1] :
          tree[index].id_R = tree_edge[i + 1];
      }

      tree[index].bl_L < 0 ?
        tree[index].bl_L = edge_length[el_index] :
        tree[index].bl_R = edge_length[el_index];
    }

    tree[0].update_path(0.0, 0.0);
  }

  std::vector<std::vector<double> > calculate_inv_path_length() {
    std::vector<std::vector<double> > res;
    for (const auto& i : tree) {
      auto inv_b_lengths = i.get_inv_blengths();
      if (!inv_b_lengths.empty()) {
        for (const auto& j : inv_b_lengths) {
          res.push_back(j);
        }
      }
    }
    return res;
  }

 private:
  std::vector< path_node > tree;
  int tree_size = 0;
};

}  // namespace mpd_tree
