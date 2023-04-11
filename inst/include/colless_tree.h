#pragma once

#include "phylotree.h"
#include "binom.h"

namespace colless {

struct colless_node {
  colless_node* daughter1 = nullptr;
  colless_node* daughter2 = nullptr;
  size_t L;
  size_t R;

  colless_node() {
    L = R = 0;
  }

  void set_both_extant() {
    L = R = 1;
  }

  void set_one_extant(colless_node& d1) {
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

struct colless_tree : public phylo_tree<colless_node> {
public:
  int calc_colless() {
    tree[0].update_num_tips();
    int s = 0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l - r < 0 ? s -= l - r : s+= l - r;
    }
    return s;
  }


  double calc_eWcolless() {
    tree[0].update_num_tips();
    double s = 0;
    for (const auto& i : tree) {
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
    for (const auto& i : tree) {
      if (i.L != i.R) s++;
    }
    return s * 1.0 / tree.size();
  }

  double calc_stairs2() {
    tree[0].update_num_tips();
    double s = 0;
    for (const auto& i : tree) {
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
        auto n1 = l;
        if (r > l) {n1 = r;}
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
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l != r ? s++ : 0;
    }
    return s;
  }

  double calc_j_one() {
    tree[0].update_num_tips();
    double s = 0.0;
    double norm_j = 0.0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      int n_j = l + r;
      norm_j += n_j;
      s += -l * std::log(1.0 * l / n_j) - r * std::log(1.0 * r / n_j);
    }
    s *= 1.0 / (norm_j * std::log(2));
    return s;
  }

  double calc_rquartet() {
    tree[0].update_num_tips();
    double s = 0.0;
    for (const auto& i : tree) {
      auto l = binom_coeff_2(i.L);  // choose_2
      auto r = binom_coeff_2(i.R);  // choose_2
      s += l * r;
    }
    return s;
  }


  double correct_pda(double Ic, size_t num_tips) {
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

  double correct_rquartet_yule(double stat, size_t num_tips) {
    auto expected = binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }

  double correct_rquartet_pda(double stat, size_t num_tips) {
    auto expected = 3.0 / 5.0 * binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }
};

}   // namespace colless
