#ifndef depth_h
#define depth_h

namespace depth {

using ltable = std::vector< std::array<double, 4>>;

class depth_stat_ltab {
public:
  depth_stat_ltab(const ltable& ltab_in) : ltable_(ltab_in) {
  }

  size_t find_parent(const ltable& ltable_,
                     int focal_id,
                     int start_index) {

    for(int i = start_index; i >= 0; i--) {
      if (static_cast<int>(ltable_[i][2]) == focal_id) {
        return i;
      }
    }

    if (start_index != ltable_.size()) {
      return find_parent(ltable_, focal_id, ltable_.size());
    } else {
      return -1; // trigger access violation --> update to throw
    }
  }

  auto collect_depths() {
    std::vector< int > s_values(ltable_.size(), 0);
    s_values[0] = 1;
    s_values[1] = 1;

    for (size_t i = 2; i < ltable_.size(); ++i) {
      int parent_index = abs(static_cast<int>(ltable_[i][1])) - 1;
      s_values[parent_index]++;
      s_values[i] = s_values[parent_index];
    }
    return s_values;
  }

  size_t calc_max_depth() {
    auto s_values = collect_depths();
    return(*std::max_element(s_values.begin(), s_values.end()));
  }

  double calc_b2() {
    auto depths = collect_depths();
    double s = 0.0;
    for (const auto& i : depths) {
      s += i / std::pow(2, i);
    }
    return s;
  }

  double calc_b1() {
    std::vector< int > depth_tracker(ltable_.size(), 1);
    double b1 = 0.0;
    std::vector<int> depths;
    for (int focal_index = ltable_.size() - 1; focal_index > 1; --focal_index) {
      auto parent = std::abs(ltable_[focal_index][1]) - 1; // parent index is in R form.
      auto max_dist = std::max(depth_tracker[focal_index],
                               depth_tracker[parent]);
      depth_tracker[parent] = 1 + max_dist;
      b1 += 1.0 / max_dist;
    }

    return b1;
  }


  double calc_var_leaf_depth() {
    std::vector< int > depths = collect_depths();

    double inv_tree_size = 1.0 / depths.size();
    double average_depth = std::accumulate(depths.begin(), depths.end(), 0.0) * inv_tree_size;

    double var_depth = 0.0;
    for (const auto& i : depths) {
      var_depth += (i - average_depth) * (i - average_depth);
    }
    var_depth *= inv_tree_size;
    return var_depth;
  }

  std::vector< int > collect_widths() {
    std::vector< int > current_depths(ltable_.size() + 1, 0);
    for (int i = 1; i < ltable_.size(); ++i) {
      int parent  = std::abs(ltable_[i][1]);
      int self_id = std::abs(ltable_[i][2]);
      current_depths.push_back(current_depths[parent]);
      current_depths[parent]++;
      current_depths[self_id] = current_depths[parent];
    }

    // make histogram
    std::vector<int> counts(ltable_.size(), 0);
    for (const auto& i : current_depths) {
      counts[i]++;
    }
    return counts;
  }

  size_t calc_max_width() {
    auto widths = collect_widths();
    return(*std::max_element(widths.begin(), widths.end()));
  }

  size_t max_del_width() {
    auto widths = collect_widths();
    std::vector<int> dW(widths.size() - 1);
    for (size_t i = 1; i < widths.size(); ++i) {
      dW[i - 1] = widths[i] - widths[i - 1];
    }
    return(*std::max_element(dW.begin(), dW.end()));
  }



private:
  const ltable ltable_;
};



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
