#ifndef colless_h
#define colless_h

#include <vector>
#include <array>
#include <numeric> // std::accumulate


using ltable = std::vector< std::array<double, 4>>;

class colless_stat_ltable {
public:
  colless_stat_ltable(const ltable& l_in) : ltable_(l_in) {
    extant_tips = std::vector<int>(l_in.size(), 1);
    num_tips = get_num_tips();
  }

  size_t calc_colless() {
    size_t colless_stat = 0;
    while(true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {// we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      colless_stat += std::abs(L - R);
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (ltable_.size() == 1) break;
    }
    return colless_stat;
  }

  size_t count_pitchforks() {
    size_t num_pitchforks = 0;
    while(true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {// we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      if (L + R == 3) num_pitchforks++;

      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (ltable_.size() == 1) break;
    }
    return num_pitchforks;
  }

  double count_stairs() {
    size_t num_s = 0;
    size_t N = ltable_.size();
    while(true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {// we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      if (L != R) num_s++;

      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (ltable_.size() == 1) break;
    }

    return num_s * 1.0 / (N - 1);
  }

  size_t count_IL() {
    size_t num_IL = 0;
    while(true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {// we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      if ((L == 1 && R > 1) ||  (L > 1 && R == 1)) {
        num_IL++;
      }

      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (ltable_.size() == 1) break;
    }
    return num_IL;
  }

  double correct_pda(double Ic) {
   double denom = powf(num_tips, 1.5f);
    return 1.0 * Ic / denom;
  }

  double correct_yule(double Ic) {
    static const double g = 0.577215664901532;
    auto output = (Ic - num_tips * log(num_tips) - num_tips * (g - 1 - log(2))) / num_tips;
    return output;
  }

private:

  size_t get_min_index() {
    auto min_val = std::min_element(ltable_.begin(), ltable_.end(),
                                    [&](const auto& a, const auto& b) {
                                      return a[0] < b[0];
                                    });
    return std::distance(ltable_.begin(), min_val);
  }

  int index_of_parent(int parent) {
    int index = 0;
    bool found = false;
    for (; index < ltable_.size(); ++index) {
      if (ltable_[index][2] == parent) {
        found = true;
        break;
      }
    }
    if (!found) index = -1;
    return index;
  }

  void remove_from_dataset(size_t index) {
    std::swap(extant_tips[index], extant_tips.back());
    extant_tips.pop_back();
    std::swap(ltable_[index], ltable_.back());
    ltable_.pop_back();
  }

  size_t get_num_tips() {
   return ltable_.size();
  }
  ltable ltable_;
  std::vector< int > extant_tips;
  size_t num_tips;
};



namespace colless_tree {


struct node {
  node* daughter1 = nullptr;
  node* daughter2 = nullptr;
  size_t L;
  size_t R;

  node() {
    L = R = 0;
  }

  void set_both_internal(node& d1, node& d2){
    daughter1 = &d1;
    daughter2 = &d2;
  }

  void set_both_extant() {
    L = R = 1;
  }

  void set_one_extant(node& d1) {
    daughter1 = &d1;
    R = 1;
  }

  size_t update_num_tips() {

    if (daughter1 && !daughter2) {
      L = daughter1->update_num_tips();

    }
    if (daughter1 && daughter2) {
      L =  daughter1->update_num_tips();
      R =  daughter2->update_num_tips();
    }

    return L + R;
  }
};


class phylo_tree {
public:

  phylo_tree(const std::vector< long >& tree_edge) {

    int root_no = static_cast<int>(tree_edge.front());
    tree.resize(tree_edge.size() / 2 - root_no + 2);

    for (size_t i = 0; i < tree_edge.size(); i += 2 ) {

      int index    = static_cast<int>(tree_edge[i]) - root_no;
      int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

      if (d1_index < 0) {
        tree[index].R == 0 ? tree[index].R = 1 : tree[index].L = 1;
      } else {
        !tree[index].daughter1 ? tree[index].daughter1 = &tree[d1_index] : tree[index].daughter2 = &tree[d1_index];
      }
    }
  }

  int calc_colless() {
    tree[0].update_num_tips();
    int s = 0;
    for(const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l - r < 0 ? s -= l - r : s+= l - r;
    }
    return s;
  }

  double calc_stairs() {
    tree[0].update_num_tips();
    int s = 0;
    for(const auto& i : tree) {
      if (i.L != i.R) s++;
    }
    return s * 1.0 / tree.size();
  }

  double correct_pda(double Ic, size_t num_tips) {
    double denom = powf(num_tips, 1.5f);
    return 1.0 * Ic / denom;
  }

  double correct_yule(double Ic, size_t num_tips) {
    static const double g = 0.577215664901532;
    auto output = (Ic - num_tips * log(num_tips) - num_tips * (g - 1 - log(2))) / num_tips;
    return output;
  }

private:
  std::vector< node > tree;
};

}





#endif
