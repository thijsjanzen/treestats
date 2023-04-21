#pragma once

#include <vector>

template <typename NODE>
class phylo_tree {
public:
  explicit phylo_tree(const std::vector< int >& tree_edge) {
    // this holds always:
    int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
    tree.resize(tree_edge.size() / 2);

    for (size_t i = 0; i < tree_edge.size(); i += 2) {
      int index    = static_cast<int>(tree_edge[i])     - root_no;
      int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

      if (d1_index > 0) {
        !tree[index].daughter1 ?
        tree[index].daughter1 = &tree[d1_index] :
        tree[index].daughter2 = &tree[d1_index];
      }
    }
  }

  std::vector< NODE > tree;
};

