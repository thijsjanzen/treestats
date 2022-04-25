#ifndef depth_h
#define depth_h

namespace depth {


struct node {
  node* daughter1 = nullptr;
  node* daughter2 = nullptr;
  size_t depth;

  node() {
    depth = 0;
  }

  size_t get_depth() {

    if (!daughter1 && !daughter2) {
      depth = 1;
    } else {
      if (daughter1 && !daughter2) {
        depth = 1 + daughter1->get_depth();
      } else {
        auto d1 = daughter1->get_depth();
        auto d2 = daughter2->get_depth();
        depth = 1 + std::max(d1, d2);
      }
    }
    return depth;
  }
};

class phylo_tree {
public:

  phylo_tree(const std::vector< long >& tree_edge) {

    int root_no = 2 + static_cast<int>(0.25 * tree_edge.size()); // this holds always.
    tree.resize(tree_edge.size() / 2);

    for (size_t i = 0; i < tree_edge.size(); i += 2 ) {

      int index    = static_cast<int>(tree_edge[i]) - root_no;
      int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

      if (d1_index > 0) {
        !tree[index].daughter1 ? tree[index].daughter1 = &tree[d1_index] : tree[index].daughter2 = &tree[d1_index];
      }
    }
  }

  int max_depth() {
    tree[0].get_depth();
    size_t md = 0;
    for (const auto& i : tree) {
      md = std::max(i.depth, md);
    }
    return md;
  }

private:
  std::vector< node > tree;
};

class depth_tracker {
public :
  depth_tracker(const std::vector< long >& tree_edge) {

    int root_no = 2 + static_cast<int>(0.25 * tree_edge.size()); // this holds always.
    tree.resize(tree_edge.size());

    for (size_t i = 0; i < tree_edge.size(); i += 2 ) {

      int index    = static_cast<int>(tree_edge[i]);
      int d1_index = static_cast<int>(tree_edge[i + 1]);

      node new_node;
      if (d1_index > root_no) { // pointing towards other node
        !tree[index].daughter1 ? tree[index].daughter1 = &tree[d1_index] : tree[index].daughter2 = &tree[d1_index];
      }
    }
  }

private:
  std::vector< node > tree;
};
}

namespace width {
struct node {
  node* daughter1 = nullptr;
  node* daughter2 = nullptr;
  int depth;
  int max_dist_to_tips;
  size_t L;
  size_t R;

  node() {
    depth = 0;
    max_dist_to_tips = 0;
  }

  int calculate_max_dist_to_tips() {
    if (!daughter1 && !daughter2) {
      max_dist_to_tips = 0;
    } else {
      if (daughter1 && !daughter2) {
        max_dist_to_tips = 1 + daughter1->calculate_max_dist_to_tips();
      } else {
        auto d1 = 1 + daughter1->calculate_max_dist_to_tips();
        auto d2 = 1 + daughter2->calculate_max_dist_to_tips();
        max_dist_to_tips = std::max(d1, d2);
      }
    }
    return max_dist_to_tips;
  }

  void set_depth(size_t parent_depth) {
    depth = 1 + parent_depth;
    if (!daughter1 && !daughter2) {

    } else {
      if (daughter1 && !daughter2) {
        daughter1->set_depth(depth);
      } else {
        daughter1->set_depth(depth);
        daughter2->set_depth(depth);
      }
    }
    return;
  }

  size_t update_l_r() {

    if (!daughter1 && !daughter2) {
      L = depth;
      R = depth;
    }

    if (daughter1 && !daughter2) {
      L = daughter1->update_l_r();
      R = depth;
    }
    if (daughter1 && daughter2) {
      L =  daughter1->update_l_r();
      R =  daughter2->update_l_r();
    }

    return L + R;
  }
};

class depth_tracker {
public :
  depth_tracker(const std::vector< long >& tree_edge) {

    root_no = 2 + static_cast<int>(0.25 * tree_edge.size()); // this holds always.
    auto max_num = *std::max_element(tree_edge.begin(), tree_edge.end());
    tree.resize(max_num + 1);

    for (size_t i = 0; i < tree_edge.size(); i += 2 ) {
      int index    = static_cast<int>(tree_edge[i]);
      int d1_index = static_cast<int>(tree_edge[i + 1]);
      !tree[index].daughter1 ? tree[index].daughter1 = &tree[d1_index] : tree[index].daughter2 = &tree[d1_index];
    }
  }

  int calc_max_width() {
    update_depth();
    std::vector<int> depths(tree.size(), 0);
    for (auto i = tree.begin() + 1; i < tree.end(); ++i) {
      depths[ (*i).depth ] ++;
    }
    return *std::max_element(depths.begin(), depths.end());
  }

  int calc_max_del_width() {
    update_depth();
    std::vector<int> depths(tree.size(), 0);
    for (auto i = tree.begin() + 1; i < tree.end(); ++i) {
      depths[ (*i).depth ] ++;
    }
    std::vector<int> dW(depths.size() - 1);
    for (size_t i = 1; i < depths.size(); ++i) {
      dW[i - 1] = depths[i] - depths[i - 1];
    }
    return(*std::max_element(dW.begin(), dW.end()));
  }

  double calc_b1() {
    tree[root_no].calculate_max_dist_to_tips();
    double b1 = 0.0;
    for (size_t i = root_no + 1; i < tree.size(); ++i) {
        b1 += 1.0 / tree[i].max_dist_to_tips;
    }
    return b1;
  }

  double calc_b2() {
    update_depth();
    double s = 0.0;
    for (size_t i = 1; i < root_no; ++i) {
      // we are only interested in tip depths
      s += tree[i].depth / std::pow(2, tree[i].depth);
    }
    return s;
  }

  double var_leaf_depth() {
    update_depth();
    double average_depth = 0;
    for (size_t i = 1; i < root_no; ++i) {
      average_depth += tree[i].depth;
    }
    average_depth *= 1.0 / (root_no - 1);
    double var_depth = 0.0;
    for (size_t i = 1; i < root_no; ++i) {
      var_depth += (tree[i].depth - average_depth) * (tree[i].depth - average_depth);
    }
    var_depth *= 1.0 / (root_no - 1);
    return var_depth;
  }

  int calc_sym_nodes() {
    update_depth();
    tree[root_no].update_l_r();
    int num_sym_nodes = 0;
    for (size_t i = root_no; i < tree.size(); ++i) {
      if (tree[i].L == tree[i].R) {
        num_sym_nodes++;
      }
    }

    return tree.size() - root_no - num_sym_nodes;
  }


private:
  void update_depth() {
    tree[root_no].set_depth(-1);
  }
  std::vector< node > tree;
  int root_no;
};



}








#endif
