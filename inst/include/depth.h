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

  size_t calc_max_depth() {
    std::vector< int > s_values(ltable_.size(), 0);
    s_values[0] = 1;
    s_values[1] = 1;

    // ltable:
    // 0 = branching time // not used here
    // 1 = parent
    // 2 = id
    // 3 = extinct time // not used here
    for (size_t i = 2; i < ltable_.size(); ++i) {
      int parent_index = abs(static_cast<int>(ltable_[i][1])) - 1;
      s_values[parent_index]++;
      s_values[i] = s_values[parent_index];
    }
    // verified with R for correct values
    // compared with apTreeShape results, based on ltable
    // 24-09-2021
    return(*std::max_element(s_values.begin(), s_values.end()));
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

}

#endif
