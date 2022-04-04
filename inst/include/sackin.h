#ifndef sackin_h
#define sackin_h

#include <vector>
#include <algorithm>
#include <array>

using ltable = std::vector< std::array<double, 4>>;

class sackin_stat_ltab {
public:
  sackin_stat_ltab(const ltable& ltab_in) : ltable_(ltab_in) {
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

  size_t calc_sackin() {
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
    return(std::accumulate(s_values.begin(), s_values.end(), 0));
  }

  double calc_blum() {
    std::vector< int > s_values(ltable_.size(), 1);
    //s_values[0] = 1;
    //s_values[1] = 1;

    // ltable:
    // 0 = branching time // not used here
    // 1 = parent
    // 2 = id
    // 3 = extinct time // not used here
    for (size_t i = ltable_.size() - 1; i > 0; i--) {
      int parent_index = abs(static_cast<int>(ltable_[i][1])) - 1;
      s_values[parent_index] += s_values[i];
      s_values[i] = s_values[parent_index];
    }

    double s = 0.0;
    for (size_t i = 1; i < s_values.size(); ++i) {
      if (s_values[i] != 0.0) {
        s += log(1.0 * s_values[i] - 1.0);
      }
    }
    return s;
  }

  double correct_pda(double Is) {
    size_t n = ltable_.size();
    double denom = powf(n, 1.5f);
    return 1.0 * Is / denom;
  }

  double correct_yule(double Is) {
    double sum_count = 0.0;
    size_t n = ltable_.size();
    for (size_t j = 2; j <= n; ++j) {
      sum_count += 1.0 / j;
    }
    return 1.0 * (Is - 2.0 * n * sum_count) / n;
  }

private:
  const ltable ltable_;
};


struct node {
  node* daughter1 = nullptr;
  node* daughter2 = nullptr;
  size_t num_extant_tips;

  node() {
    num_extant_tips = 0;
  }

  size_t get_acc_num_tips() {

    if (!daughter1 && !daughter2) {
      num_extant_tips = 2;
    } else {
      if (daughter1 && !daughter2) {
        num_extant_tips = 1 + daughter1->get_acc_num_tips();
      } else {
        num_extant_tips = daughter1->get_acc_num_tips() + daughter2->get_acc_num_tips();
      }
    }
    return num_extant_tips;
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

  int calc_sackin() {
    tree[0].get_acc_num_tips();
    int s = 0;
    for(const auto& i : tree) {
      s += i.num_extant_tips;
    }
    return s;
  }

  size_t count_pitchforks() {
    tree[0].get_acc_num_tips();
    size_t num_pitchforks = 0;
    for(const auto& i : tree) {
      if (i.num_extant_tips == 3) {
        num_pitchforks++;
      }
    }
    return num_pitchforks;
  }

  size_t count_cherries() {
    tree[0].get_acc_num_tips();
    size_t num_cherries = 0;
    for(const auto& i : tree) {
      if (i.num_extant_tips == 2) {
        num_cherries++;
      }
    }
    return num_cherries;
  }

  double correct_pda(size_t n,
                     double Is) {
    double denom = powf(n, 1.5f);
    return 1.0 * Is / denom;
  }

  double correct_yule(size_t n,
                      double Is) {
    double sum_count = 0.0;
    for (size_t j = 2; j <= n; ++j) {
      sum_count += 1.0 / j;
    }
    return 1.0 * (Is - 2.0 * n * sum_count) / n;
  }

  double calc_blum() {
    tree[0].get_acc_num_tips();
    double s = 0;
    for(const auto& i : tree) {
      if (i.num_extant_tips > 1) {
        s += log(1.0 * i.num_extant_tips - 1);
      }
    }
    return s;
  }

private:
  std::vector< node > tree;
};


#endif
