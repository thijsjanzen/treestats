#ifndef colless_h
#define colless_h

#include <vector>
#include <array>
#include <numeric> // std::accumulate

#include <iostream>

enum analysis_type {colless, ewColless, rogers, pitchforks, stairs2,
                      IL};

using ltable = std::vector< std::array<double, 4>>;

class colless_stat_ltable {
public:
  colless_stat_ltable(const ltable& l_in) : ltable_(l_in) {
    extant_tips = std::vector<int>(l_in.size(), 1);
    num_tips = get_num_tips();
  }

  template<typename ANAL_TYPE>
  constexpr double update_stat(int L, int R, double stat) {

    if constexpr (ANAL_TYPE == analysis_type::colless) {
      return stat + std::abs(L + R);
    }

    if constexpr (ANAL_TYPE == analysis_type::ewColless) {
      if (L + R > 2) {
        stat += 1.0 * std::abs(L - R) / (L + R - 2);
      }
      return stat;
    }

    if constexpr (ANAL_TYPE == analysis_type::rogers) {
      if (L != R) stat++;
      return stat;
    }

    if constexpr (ANAL_TYPE == analysis_type::pitchforks) {
      if (L + R == 3) stat++;
      return stat;
    }

    if constexpr (ANAL_TYPE == analysis_type::stairs2) {
      int min_l_r, max_l_r;
      if (L < R) {
        min_l_r = L; max_l_r = R;
      } else {
        min_l_r = R; max_l_r = L;
      }
      stat += 1.0 * min_l_r / max_l_r;
      return stat;
    }

    if constexpr (ANAL_TYPE == analysis_type::IL) {
      if ((L == 1 && R > 1) ||  (L > 1 && R == 1)) {
        stat++;
      }
      return stat;
    }

  }


  template<typename ANAL_TYPE>
  double collect_stat() {
    double stat = 0.0;
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
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      stat = update_stat<ANAL_TYPE>(L, R, stat);

      if (ltable_.size() == 1) break;
    }
    return stat;
  }

  size_t calc_colless() {
    size_t colless_stat = collect_stat<analysis_type::colless>();
    return colless_stat;
  }

  double calc_ew_colless() {
    int N = ltable_.size();
    if (N <= 2) return 0;
    double ew_colless_stat = collect_stat<analysis_type::ewColless>();
    return ew_colless_stat * 1.0 / (N - 2);
  }


  size_t calc_rogers() {
    size_t rogers_stat = collect_stat<analysis_type::rogers>();
    return rogers_stat;
  }

  size_t count_pitchforks() {
    size_t num_pitchforks = collect_stat<analysis_type::pitchforks>();
    return num_pitchforks;
  }

  double count_stairs() {
    size_t num_s = collect_stat<analysis_type::rogers>();
    size_t N = ltable_.size();
    return num_s * 1.0 / (N - 1);
  }

  double count_stairs2() {
    double num_s = collect_stat<analysis_type::stairs2>();
    size_t N = ltable_.size();
    return num_s * 1.0 / (N - 1);
  }

  double count_IL() {
    return collect_stat<analysis_type::IL>();
  }

  std::vector<double> collect_I() {
    std::vector<double> i_vals;
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

      int L_R = L + R;
      if (L_R > 3) {
        double avg_n = std::ceil(L_R * 0.5); // N / 2 + N % 2 (see Fusco 1995).
        double I_val =  1.0 * (std::max(L, R) - avg_n) / ((L_R - 1) - avg_n);
        if (L_R % 2 == 0) {
          I_val *= 1.0 * (L_R - 1) / L_R;
        }
        i_vals.push_back(I_val);
      }

      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (ltable_.size() == 1) break;
    }
    return i_vals;
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

    int root_no = 2 + static_cast<int>(0.25 * tree_edge.size()); // this holds always.

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


  double calc_eWcolless() {
    tree[0].update_num_tips();
    double s = 0;
    for(const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      double l_r = l + r;
      if (l_r > 2) {
        s += std::abs(l - r) * 1.0 / (l_r - 2);
      }
    }
    s *= 1.0 / (tree.size() - 1);
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

  double calc_stairs2() {
    tree[0].update_num_tips();
    double s = 0;
    for(const auto& i : tree) {
      int min_l_r, max_l_r;
      if (i.L < i.R) {
        min_l_r = i.L; max_l_r = i.R;
      } else {
        min_l_r = i.R; max_l_r = i.L;
      }
      s += 1.0 * min_l_r / max_l_r;
    }
    return s * 1.0 / tree.size();
  }

  std::vector<double> collect_I() {
    tree[0].update_num_tips();
    std::vector<double> i_vals;
    for (size_t i = 0; i < tree.size(); ++i) {
      int l = tree[i].L;
      int r = tree[i].R;
      int nv = l + r;
      if (nv > 3) {
        double avg_n = std::ceil(nv * 0.5);
        auto n1 = l; if (r > l) n1 = r;
        double I_val =  1.0 * (n1 - avg_n) / ((nv - 1) - avg_n);

        if (nv % 2 == 0) {
          I_val *= 1.0 * (nv - 1) / nv;
        }
        i_vals.push_back(I_val);
      }
    }
    return i_vals;
  }

  int calc_rogers() {
    tree[0].update_num_tips();
    int s = 0;
    for(const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l != r ? s++ : 0;
    }
    return s;
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
