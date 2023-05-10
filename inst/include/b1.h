
#pragma once

#include "phylotree.h"

namespace b1_tree {
  class b1_tree {
    struct node_t {
      node_t* daughter1 = nullptr;
      node_t* daughter2 = nullptr;
      int max_dist_to_tips = 1;

      void update_max_dist_to_tips() {
        if (daughter1 && !daughter2) {
            max_dist_to_tips = 1 + daughter1->max_dist_to_tips;
          } else if (daughter1 && daughter2) {
            auto d1 = 1 + daughter1->max_dist_to_tips;
            auto d2 = 1 + daughter2->max_dist_to_tips;
            max_dist_to_tips = std::max(d1, d2);
        }
        return;
      }
    };

    phylo_tree_t<node_t> tree;

  public:
    explicit b1_tree(const std::vector< int >& tree_edge)
      : tree(make_phylo_tree<node_t>(tree_edge)) {
    }

    double calc_b1() {
      double b1 = 0.0;
      for (size_t i = tree.size() - 1; i >= 1; --i) {
        tree[i].update_max_dist_to_tips();
        b1 += 1.0 / tree[i].max_dist_to_tips;
      }
      return b1;
    }
  };

} // namespace b1_tree
