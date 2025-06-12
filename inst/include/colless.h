// Copyright 2022 - 2025 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
#pragma once

#include <vector>
#include <array>
#include <numeric>  // std::accumulate
#include <utility>  // swap
#include <algorithm>  // std::min_element

#include "binom.h"       // NOLINT [build/include_subdir]
#include "phylotree.h"   // NOLINT [build/include_subdir]

double calc_colless(int L, int R) {
  return std::abs(L - R);
}

double calc_colless_quad(int L, int R) {
  return (L - R) * (L - R);
}

double root_imbal(int L, int R) {
  return L + R;
}

double calc_ew_colless(int L, int R) {
  double answ = 0.0;
  if (L + R > 2) {
    answ = std::abs(L - R) * 1.0 / (L + R - 2);
  }
  return answ;
}

double calc_stairs(int L, int R) {
  if (L != R) return 1;
  return 0;
}

double calc_stairs2(int L, int R) {
  if (L < R) {
    return 1.0 * L / R;
  } else {
    return 1.0 * R / L;
  }
}

double calc_rogers(int L, int R) {
  if (L != R) return 1;
  return 0;
}

double calc_j_one(int L, int R) {
  int n_j = L + R;
  return -L * std::log(1.0 * L / n_j) - R * std::log(1.0 * R / n_j);
}

double calc_rquartet(int L, int R) {
  return binom_coeff_2(L) * binom_coeff_2(R);
}

double calc_I(int L, int R) {
  double I_val = 0.0;
  int nv = L + R;
  if (nv > 3) {
    double avg_n = std::ceil(nv * 0.5);

    auto n1 = R > L ? R : L;

    I_val =  1.0 * (n1 - avg_n) / ((nv - 1) - avg_n);
    if (nv % 2 == 0) {
      I_val *= 1.0 * (nv - 1) / nv;
    }
  }
  return I_val;
}

double calc_ILnumber(int L, int R) {
  if ((L == 1 && R > 1) ||  (L > 1 && R == 1)) {
    return 1.0;
  }
  return 0.0;
}

double calc_pitchforks(int L, int R) {
  if (L + R == 3) return 1.0;
  return 0.0;
}

namespace colless_tree {

class colless_tree {
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    size_t L = 1;    // cached value, valid after 'update_num_tips()'
    size_t R = 1;    // cached value, valid after 'update_num_tips()'

    double update_node(double (*calculation)(int, int)) {
      if (daughter1 && !daughter2) {
        L = daughter1->L + daughter1->R;
      } else if (daughter1 && daughter2) {
        L = daughter1->L + daughter1->R;
        R = daughter2->L + daughter2->R;
      }

      return calculation(L, R);
    }
  };

  phylo_tree_t<node_t> tree;

 public:
  explicit colless_tree(const std::vector< int >& tree_edge)
    : tree(make_phylo_tree<node_t, false>(tree_edge)) {
  }

  double calc_stat(double (*calculation)(int, int)) {
    double s = 0.0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      s += (*i).update_node(calculation);
    }
    return s;
  }

  double calc_j_one(double (*calculation)(int, int)) {
    double s = 0.0;
    double norm_j = 0.0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      s += (*i).update_node(calculation);
      norm_j += (*i).L + (*i).R;
    }
    s *= 1.0 / (norm_j * std::log(2));
    return s;
  }

  double collect_I(double (*calculation)(int, int)) {
    double s = 0.0;
    int cnt = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      double val = (*i).update_node(calculation);
      if ((*i).L + (*i).R > 3) {
        cnt++;
        s += val;
      }
    }
    return s * 1.0 / cnt;
  }

  double calc_root_imbal() {
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      (*i).update_node(&root_imbal);
    }
    auto root = tree.begin();
    auto n1 = root->L;
    auto n2 = root->R;
    double answ = 1.0 * n1 / (n1 + n2);
    if (answ < 0.5) answ = 1.0 - answ;
    return answ;
  }

  size_t calc_double_cherries() {
    size_t num = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      (*i).update_node(&root_imbal);

      if ((*i).L == 2 && (*i).R == 2) {
        num++;
      }
    }
    return num;
  }

  size_t calc_four_prong() {
    size_t num = 0;
    for (auto i = tree.rbegin(); i != tree.rend(); ++i) {
      (*i).update_node(&root_imbal);

      if ((*i).L == 3 && (*i).R == 1) {
        num++;
      } else if ((*i).L == 1 && (*i).R == 3) {
        num++;
      }
    }
    return num;
  }

  int size() {
    return tree.size();
  }
};

}    // end namespace colless_tree

using ltable = std::vector< std::array<double, 4>>;

class colless_stat_ltable {
 public:
  explicit colless_stat_ltable(const ltable& l_in) : ltable_(l_in) {
    extant_tips = std::vector<int>(l_in.size(), 1);
    num_tips = get_num_tips();
  }

  double collect_stat(double (*calculation)(int, int)) {
    double stat = 0.0;
    while (true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {  // we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      stat += calculation(L, R);

      if (ltable_.size() == 1) break;
    }
    return stat;
  }

  size_t colless() {
    return collect_stat(&calc_colless);
  }

  size_t colless_quad() {
    return collect_stat(&calc_colless_quad);
  }

  double ew_colless() {
    int N = ltable_.size();
    if (N <= 2) return 0;
    double ew_colless_stat = collect_stat(&calc_ew_colless);
    return ew_colless_stat * 1.0 / (N - 2);
  }

  size_t rogers() {
    size_t rogers_stat = collect_stat(&calc_rogers);
    return rogers_stat;
  }

  size_t count_pitchforks() {
    return collect_stat(&calc_pitchforks);
  }

  double count_stairs() {
    size_t N = ltable_.size();
    size_t num_s = collect_stat(&calc_stairs);
    return num_s * 1.0 / (N - 1);
  }

  double count_stairs2() {
    size_t N = ltable_.size();
    double num_s = collect_stat(&calc_stairs2);
    return num_s * 1.0 / (N - 1);
  }

  double count_IL() {
    return collect_stat(&calc_ILnumber);
  }

  double count_rquartet() {
    return collect_stat(&calc_rquartet);
  }

  double collect_I() {
    double s_i_val = 0.0;
    int cnt = 0;
    while (true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {  // we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];

      int L_R = L + R;
      s_i_val += calc_I(L, R);
      if (L_R > 3) cnt++;

      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (ltable_.size() == 1) break;
    }
    return s_i_val * 1.0 / cnt;
  }

  double collect_j_one() {
    double stat = 0.0;
    double sum_nj = 0.0;
    while (true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {  // we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      double l_r = L + R;
      sum_nj += l_r;
      stat += calc_j_one(L, R);

      if (ltable_.size() == 1) break;
    }
    stat *= 1.0 / (sum_nj * std::log(2));
    return stat;
  }

  double calc_double_cherries() {
    double stat = 0.0;
    while (true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {  // we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (L == 2 && R == 2) stat++;

      if (ltable_.size() == 1) break;
    }
    return stat;
  }

  double calc_four_prong() {
    double stat = 0.0;
    while (true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {  // we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (L == 3 && R == 1) stat++;
      if (L == 1 && R == 3) stat++;

      if (ltable_.size() == 1) break;
    }
    return stat;
  }

  size_t size() {
    return num_tips;
  }


 private:
  size_t get_min_index() {
    auto min_val = std::min_element(ltable_.begin(), ltable_.end(),
                                    [&](const auto& a, const auto& b) {
                                      return a[0] < b[0];
                                    });
    return std::distance(ltable_.begin(), min_val);
  }

  size_t index_of_parent(int parent) {
    size_t index = 0;
    bool found = false;
    for (; index < ltable_.size(); ++index) {
      if (ltable_[index][2] == parent) {
        found = true;
        break;
      }
    }
    if (!found) {
      throw "can't find parent\n";
    }
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


// correction functions:
namespace correction {
double correct_pda(double Ic, size_t num_tips) {  // colless corrections
    double denom = powf(num_tips, 1.5f);
    return 1.0 * Ic / denom;
}

double correct_yule(double Ic, size_t num_tips) {
    static const double g = 0.577215664901532;
    auto output = (Ic -
                   num_tips * log(num_tips) -
                   num_tips * (g - 1 - log(2))) / num_tips;
    return output;
}

}   // namespace correction
