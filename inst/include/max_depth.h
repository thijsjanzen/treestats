#pragma once
#include "phylotree.h"

namespace max_depth {
class max_depth_tree {
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    size_t depth = 0;

    void set_depth() {
      if (!daughter1 && !daughter2) {
        depth = 1;
      } else {
        if (daughter1 && !daughter2) {
          depth = 1 + daughter1->depth;
        } else {
          auto d1 = daughter1->depth;
          auto d2 = daughter2->depth;
          depth = 1 + std::max(d1, d2);
        }
      }
    }
  };

  phylo_tree_t<node_t> tree;

public:
  explicit max_depth_tree(const std::vector< int >& tree_edge)
    : tree(make_phylo_tree<node_t>(tree_edge)) {
  }

  int max_depth() {
    size_t md = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      (*i).set_depth();
      if ((*i).depth > md) md = (*i).depth;
    }
    return md;
  }
};

}  // end namespace max_depth

