#pragma once

#include <cmath>
#include <algorithm>
#include "phylotree.h"
#include "binom.h"


namespace colless {


namespace original {

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

}


namespace solution1  {

struct colless_node {
  colless_node* daughter1 = nullptr;
  colless_node* daughter2 = nullptr;
  size_t L = 1;    // cached value, valid after 'update_num_tips()'
  size_t R = 1;    // cached value, valid after 'update_num_tips()'

  size_t update_num_tips() noexcept {
    L = (daughter1) ? daughter1->update_num_tips() : 1;
    R = (daughter2) ? daughter2->update_num_tips() : 1;
    return L + R;
  }
};


struct colless_tree : public phylo_tree<colless_node> {
public:
  explicit colless_tree(const std::vector< int >& tree_edge) 
  : phylo_tree<colless_node>(tree_edge) {
  }

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


}   // namespace solution1


namespace solution2  {

// 'phylotree.h

template <typename NODE>
using phylo_tree_t = std::vector<NODE>;


// standalone 'virtual contructor'
template <typename NODE>
auto make_phylo_tree(const std::vector<int> tree_edge) {
  // this holds always:
  int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
  auto tree = phylo_tree_t<NODE>(tree_edge.size() / 2);

  for (size_t i = 0; i < tree_edge.size(); i += 2) {
    int index    = static_cast<int>(tree_edge[i])     - root_no;
    int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

    if (d1_index > 0) {
      !tree[index].daughter1 ?
      tree[index].daughter1 = &tree[d1_index] :
      tree[index].daughter2 = &tree[d1_index];
    }
  }
  return tree;
}

// 'phylotree.h ends


class colless_tree 
{
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    size_t L = 1;    // cached value, valid after 'update_num_tips()'
    size_t R = 1;    // cached value, valid after 'update_num_tips()'

    size_t update_num_tips() noexcept {
      L = (daughter1) ? daughter1->update_num_tips() : 1;
      R = (daughter2) ? daughter2->update_num_tips() : 1;
      return L + R;
    }
  };

  phylo_tree_t<node_t> tree;

public:
  explicit colless_tree(const std::vector< int >& tree_edge) 
  : tree(make_phylo_tree<node_t>(tree_edge)) {
    tree[0].update_num_tips();
  }


  int calc_colless() const {
    int s = 0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l - r < 0 ? s -= l - r : s+= l - r;
    }
    return s;
  }


  double calc_eWcolless() const {
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

  double calc_stairs() const {
    int s = 0;
    for (const auto& i : tree) {
      if (i.L != i.R) s++;
    }
    return s * 1.0 / tree.size();
  }

  double calc_stairs2() const {
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

  std::vector<double> collect_I() const {
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

  int calc_rogers() const {
    int s = 0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l != r ? s++ : 0;
    }
    return s;
  }

  double calc_j_one() const {
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

  double calc_rquartet() const {
    double s = 0.0;
    for (const auto& i : tree) {
      auto l = binom_coeff_2(i.L);  // choose_2
      auto r = binom_coeff_2(i.R);  // choose_2
      s += l * r;
    }
    return s;
  }


  double correct_pda(double Ic, size_t num_tips) const {
    double denom = powf(num_tips, 1.5f);
    return 1.0 * Ic / denom;
  }

  double correct_yule(double Ic, size_t num_tips) const {
    static const double g = 0.577215664901532;
    auto output = (Ic -
                   num_tips * log(num_tips) -
                   num_tips * (g - 1 - log(2))) / num_tips;
    return output;
  }

  double correct_rquartet_yule(double stat, size_t num_tips) const {
    auto expected = binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }

  double correct_rquartet_pda(double stat, size_t num_tips) const {
    auto expected = 3.0 / 5.0 * binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }
};


}   // namespace solution2


namespace solution3  {

// 'phylotree.h

template <typename NODE>
using phylo_tree_t = std::vector<NODE>;


// standalone 'virtual contructor'
template <typename NODE>
auto make_phylo_tree(const std::vector<int> tree_edge) {
  // this holds always:
  int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
  auto tree = phylo_tree_t<NODE>(tree_edge.size() / 2);

  for (size_t i = 0; i < tree_edge.size(); i += 2) {
    int index    = static_cast<int>(tree_edge[i])     - root_no;
    int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

    if (d1_index > 0) {
      !tree[index].daughter1 ?
      tree[index].daughter1 = &tree[d1_index] :
      tree[index].daughter2 = &tree[d1_index];
    }
  }
  return tree;
}

// 'phylotree.h ends


class colless_tree 
{
  struct node_t {
    node_t* daughter1 = nullptr;
    node_t* daughter2 = nullptr;
    size_t L = 1;    // cached value, valid after 'update_num_tips()'
    size_t R = 1;    // cached value, valid after 'update_num_tips()'
  };

  phylo_tree_t<node_t> tree;


  size_t update_num_tips(node_t* node) const noexcept {
    node->L = (node->daughter1) ? update_num_tips(node->daughter1) : 1;
    node->R = (node->daughter2) ? update_num_tips(node->daughter2) : 1;
    return node->L + node->R;
  }


public:
  explicit colless_tree(const std::vector< int >& tree_edge) 
  : tree(make_phylo_tree<node_t>(tree_edge)) {
    update_num_tips(&tree[0]);
  }


  int calc_colless() const {
    int s = 0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l - r < 0 ? s -= l - r : s+= l - r;
    }
    return s;
  }


  double calc_eWcolless() const {
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

  double calc_stairs() const {
    int s = 0;
    for (const auto& i : tree) {
      if (i.L != i.R) s++;
    }
    return s * 1.0 / tree.size();
  }

  double calc_stairs2() const {
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

  std::vector<double> collect_I() const {
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

  int calc_rogers() const {
    int s = 0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l != r ? s++ : 0;
    }
    return s;
  }

  double calc_j_one() const {
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

  double calc_rquartet() const {
    double s = 0.0;
    for (const auto& i : tree) {
      auto l = binom_coeff_2(i.L);  // choose_2
      auto r = binom_coeff_2(i.R);  // choose_2
      s += l * r;
    }
    return s;
  }


  double correct_pda(double Ic, size_t num_tips) const {
    double denom = powf(num_tips, 1.5f);
    return 1.0 * Ic / denom;
  }

  double correct_yule(double Ic, size_t num_tips) const {
    static const double g = 0.577215664901532;
    auto output = (Ic -
                   num_tips * log(num_tips) -
                   num_tips * (g - 1 - log(2))) / num_tips;
    return output;
  }

  double correct_rquartet_yule(double stat, size_t num_tips) const {
    auto expected = binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }

  double correct_rquartet_pda(double stat, size_t num_tips) const {
    auto expected = 3.0 / 5.0 * binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }
};


}   // namespace solution3


namespace solution4  {

// 'phylotree.h

template <typename T>
struct node_t {
  T* daughter1 = nullptr;
  T* daughter2 = nullptr;
};


template <typename NODE>
using tree_t = std::vector<NODE>;


template <typename NODE>
tree_t<NODE> make_tree(const std::vector<int> tree_edge) {
  // this holds always:
  int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
  auto tree = tree_t<NODE>(tree_edge.size() / 2);

  for (size_t i = 0; i < tree_edge.size(); i += 2) {
    int index    = static_cast<int>(tree_edge[i])     - root_no;
    int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

    if (d1_index > 0) {
      !tree[index].daughter1 ?
      tree[index].daughter1 = &tree[d1_index] :
      tree[index].daughter2 = &tree[d1_index];
    }
  }
  return tree;
}


template <typename NLHS, typename NRHS>
tree_t<NLHS> tree_cast(const tree_t<NRHS>& rhs) {
  auto lhs = tree_t<NLHS>(rhs.size());
  auto* lhs_first = lhs.data();
  auto* rhs_first = rhs.data();
  for (size_t i = 0; i < rhs.size(); ++i) {
    lhs[i].daughter1 = (rhs[i].daughter1) ? lhs_first + std::distance<const NRHS*>(rhs_first, rhs[i].daughter1) : nullptr;
    lhs[i].daughter2 = (rhs[i].daughter2) ? lhs_first + std::distance<const NRHS*>(rhs_first, rhs[i].daughter2) : nullptr;
  }
  return lhs;
}

// 'phylotree.h ends

// augmented node
struct colless_node_t : public node_t<colless_node_t> {
  size_t L = -1;
  size_t R = -1;
};


class colless_tree 
{
  tree_t<colless_node_t> tree;

  size_t update_num_tips(colless_node_t* node) const noexcept {
    node->L = (node->daughter1) ? update_num_tips(node->daughter1) : 1;
    node->R = (node->daughter2) ? update_num_tips(node->daughter2) : 1;
    return node->L + node->R;
  }


public:
  template <typename NRHS>
  explicit colless_tree(const tree_t<NRHS>& rhs) 
  : tree(tree_cast<colless_node_t>(rhs)) {
    update_num_tips(&tree[0]);
  }


  explicit colless_tree(const std::vector<int>& tree_edge) 
  : colless_tree(make_tree<colless_node_t>(tree_edge)) {
  }


  int calc_colless() const {
    int s = 0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l - r < 0 ? s -= l - r : s+= l - r;
    }
    return s;
  }


  double calc_eWcolless() const {
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

  double calc_stairs() const {
    int s = 0;
    for (const auto& i : tree) {
      if (i.L != i.R) s++;
    }
    return s * 1.0 / tree.size();
  }

  double calc_stairs2() const {
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

  std::vector<double> collect_I() const {
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

  int calc_rogers() const {
    int s = 0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      l != r ? s++ : 0;
    }
    return s;
  }

  double calc_j_one() const {
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

  double calc_rquartet() const {
    double s = 0.0;
    for (const auto& i : tree) {
      auto l = binom_coeff_2(i.L);  // choose_2
      auto r = binom_coeff_2(i.R);  // choose_2
      s += l * r;
    }
    return s;
  }


  double correct_pda(double Ic, size_t num_tips) const {
    double denom = powf(num_tips, 1.5f);
    return 1.0 * Ic / denom;
  }

  double correct_yule(double Ic, size_t num_tips) const {
    static const double g = 0.577215664901532;
    auto output = (Ic -
                   num_tips * log(num_tips) -
                   num_tips * (g - 1 - log(2))) / num_tips;
    return output;
  }

  double correct_rquartet_yule(double stat, size_t num_tips) const {
    auto expected = binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }

  double correct_rquartet_pda(double stat, size_t num_tips) const {
    auto expected = 3.0 / 5.0 * binom_coeff(num_tips, 4);
    return stat * 1.0 / expected;
  }
};


}   // namespace solution5

using solution4::colless_tree;

}   // namespace colless
