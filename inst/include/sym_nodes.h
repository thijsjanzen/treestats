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
#include <algorithm>
#include <vector>

#include "phylotree.h"   // NOLINT [build/include_subdir]

namespace sym_nodes {

class sym_node_tree {
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    int depth = 0;
    size_t L;
    size_t R;
    std::vector<size_t> L_vec;
    std::vector<size_t> R_vec;

    void set_depth(size_t parent_depth) {
      depth = 1 + parent_depth;

      if (daughter1 && daughter2) {
        daughter1->set_depth(depth);
        daughter2->set_depth(depth);
      }

      return;
    }

    size_t update_l_r() {
      L = daughter1 ? daughter1->update_l_r() : depth;
      R = daughter2 ? daughter2->update_l_r() : depth;

      return L + R;
    }

    std::vector<size_t> update_vecs() {
      std::vector<size_t> L_R_vec;
      if (!daughter1 && !daughter2) {
        L_vec = {L};
        R_vec = {R};
      }

      if (daughter1 && daughter2) {
        L_vec = daughter1->update_vecs();
        R_vec = daughter2->update_vecs();
        L = std::accumulate(L_vec.begin(), L_vec.end(), 0.0) + depth;
        R = std::accumulate(R_vec.begin(), R_vec.end(), 0.0) + depth;
      }

      L_R_vec = L_vec;
      L_R_vec.insert(L_R_vec.end(), R_vec.begin(), R_vec.end());

      L_R_vec.push_back(L);
      L_R_vec.push_back(R);

      return L_R_vec;
    }
  };

  void update_depth() {
    tree[root_no].set_depth(-1);
  }

  bool compare_depth_dist(std::vector<size_t>* v1,
                          std::vector<size_t>* v2) {
    if ((*v1).size() != (*v2).size()) return false;

    std::sort((*v1).begin(), (*v1).end());
    std::sort((*v2).begin(), (*v2).end());

    for (size_t i = 0; i < (*v1).size(); ++i) {
      if ((*v1)[i] != (*v2)[i]) return false;
    }
    return true;
  }
  phylo_tree_t<node_t> tree;
  int root_no;

 public:
  explicit sym_node_tree(const std::vector< int >& tree_edge)
    : tree(make_phylo_tree<node_t, true>(tree_edge)) {
    root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
  }

  int calc_sym_nodes() {
    update_depth();
    tree[root_no].update_l_r();
    tree[root_no].update_vecs();
    int num_sym_nodes = 0;
    for (size_t i = root_no; i < tree.size(); ++i) {
      if (tree[i].L == tree[i].R) {
        if (compare_depth_dist(&tree[i].L_vec,
                               &tree[i].R_vec)) {
          num_sym_nodes++;
        }
      }
    }

    return tree.size() - root_no - num_sym_nodes;
  }
};

}   // namespace sym_nodes
