# pragma once

#include <vector>
#include <array>
#include <cmath>

namespace branch_colless {


// here, we don't use phylo_tree.h, because we also need to
// track edge lengths. Code is very similar though.
struct bc_node {
  bc_node* daughterL = nullptr;
  bc_node* daughterR = nullptr;
  double bl_R;
  double bl_L;
  double sum_bl_R; // sum branch lenghts connected to daughter R
  double sum_bl_L; // sum branch lenghts connected to daughter L

  bc_node() {
    daughterL = nullptr;
    daughterR = nullptr;
    bl_R = -1.0;
    bl_L = -1.0;
    sum_bl_R = 0.0;
    sum_bl_R = 0.0;
  }

  double update_node_info() {
    sum_bl_R = bl_R;
    sum_bl_L = bl_L;
    if (daughterL && !daughterR) {
      sum_bl_L = bl_L + daughterL->update_node_info();
    }
    if (daughterR && !daughterL) {
      sum_bl_R = bl_R + daughterR->update_node_info();
    }
    if (daughterL && daughterR) {
      sum_bl_L =  bl_L + daughterL->update_node_info();
      sum_bl_R =  bl_R + daughterR->update_node_info();
    }

    return sum_bl_L + sum_bl_R;
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

    tree[0].update_node_info();
  }

  double calculate_branch_colless() {
    double stat = 0.0;
    for (const auto& i : tree) {
      stat += std::fabs(i.sum_bl_L - i.sum_bl_R);
    }
    return stat;
  }

private:
  std::vector< bc_node > tree;
  int tree_size = 0;
};
} // end namespace branch_colless
