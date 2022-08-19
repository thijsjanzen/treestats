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

}

namespace width {
struct node {
  node* daughter1 = nullptr;
  node* daughter2 = nullptr;
  int depth;
  int max_dist_to_tips;
  size_t L;
  size_t R;
  std::vector<size_t> L_vec;
  std::vector<size_t> R_vec;

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
    if (daughter1 && !daughter2) {
      daughter1->set_depth(depth);
    }

    if (daughter2 && !daughter1) {
      daughter2->set_depth(depth);
    }

    if (daughter1 && daughter2) {
      daughter1->set_depth(depth);
      daughter2->set_depth(depth);
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

  std::vector<size_t> update_vecs() {

    std::vector<size_t> L_R_vec;

    if (!daughter1 && !daughter2) {
      L_vec = {L};
      R_vec = {R};
    }

    if (daughter1 && !daughter2) {
      L_vec = daughter1->update_vecs();
      L = std::accumulate(L_vec.begin(), L_vec.end(), 0.0) + depth;
      R_vec = {R};
    }

    if (daughter2 && !daughter1) {
      R_vec = daughter2->update_vecs();
      R = std::accumulate(R_vec.begin(), R_vec.end(), 0.0) + depth;
      L_vec = {L};
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
    double average_depth = 0.0;

    int n = root_no - 1;
    for (size_t i = 1; i < root_no; ++i) {
      average_depth += tree[i].depth;

  //    std::cerr << tree[i].depth << "\n";
    }

    average_depth *= 1.0 / (n);
  //  std::cerr << average_depth << " " << root_no - 1 << "\n";

    double var_depth = 0.0;
    for (size_t i = 1; i < root_no; ++i) {
      var_depth += (tree[i].depth - average_depth) * (tree[i].depth - average_depth);
    }
    var_depth *= 1.0 / (n - 1);
    return var_depth;
  }

  int calc_sym_nodes() {
    update_depth();
    tree[root_no].update_l_r();
    tree[root_no].update_vecs();
    int num_sym_nodes = 0;
    for (size_t i = root_no; i < tree.size(); ++i) {
      if (tree[i].L == tree[i].R) {
        if (compare_depth_dist(tree[i].L_vec, tree[i].R_vec)) {
          num_sym_nodes++;
        }
      }
    }

    return tree.size() - root_no - num_sym_nodes;
  }

private:
  void update_depth() {
    tree[root_no].set_depth(-1);
  }

  bool compare_depth_dist(std::vector<size_t>& v1,
                          std::vector<size_t>& v2) {
    if (v1.size() != v2.size()) return false;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    for (size_t i = 0; i < v1.size(); ++i) {
      if (v1[i] != v2[i]) return false;
    }
    return true;
  }
  std::vector< node > tree;
  int root_no;
};



}

#endif
