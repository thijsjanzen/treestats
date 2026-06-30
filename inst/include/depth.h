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

#include <algorithm>   // std::sort
#include <memory>
#include <vector>
#include "util.h"      // NOLINT [build/include_subdir]
#include "binom.h"     // NOLINT [build/include_subdir]


namespace depth {

struct node_binary {
  node_binary* daughter1 = nullptr;
  node_binary* daughter2 = nullptr;

  int depth;
  int dist_to_tips;     // for b1 statistic
  int num_extant_tips;  // for sackin related statistics

  node_binary() {
    depth = 0;
    dist_to_tips = 0;
    num_extant_tips = 0;
    daughter1 = nullptr;
    daughter2 = nullptr;
  }

  void add_daughter(node_binary* d) {
    if (!daughter1) {
      daughter1 = d;
    } else {
      daughter2 = d;
    }
  }

  void set_depth(int parent_depth) {
    depth = 1 + parent_depth;

    if (daughter1 && daughter2) {
      daughter1->set_depth(depth);
      daughter2->set_depth(depth);
    }

    return;
  }

  int set_dist_to_tips() {
    if (daughter1 && !daughter2) {
      dist_to_tips = 1 + daughter1->set_dist_to_tips();
    } else if (daughter1 && daughter2) {
      auto d1 = 1 + daughter1->set_dist_to_tips();
      auto d2 = 1 + daughter2->set_dist_to_tips();
      dist_to_tips = std::max(d1, d2);
    }

    return dist_to_tips;
  }

  int get_acc_num_tips() {
    if (!daughter1 && !daughter2) {
      num_extant_tips = 1;
    } else {
      if (daughter1 && !daughter2) {
        num_extant_tips = 1 + daughter1->get_acc_num_tips();
      } else {
        num_extant_tips = daughter1->get_acc_num_tips() +
          daughter2->get_acc_num_tips();
      }
    }
    return num_extant_tips;
  }


  bool is_double_cherry() {
    if (num_extant_tips != 4) return false;

    int L = daughter1->num_extant_tips;
    int R = daughter2->num_extant_tips;
    if (L == 2 && R == 2) return true;

    return false;
  }

  bool is_four_prong() {
    if (num_extant_tips != 4) return false;

    int L = daughter1->num_extant_tips;
    int R = daughter2->num_extant_tips;
    if (L == 1 && R == 3) return true;
    if (L == 3 && R == 1) return true;

    return false;
  }

  double calc_root_imbalance() {
    int L = daughter1->num_extant_tips;

    double answ = 1.0 * L / num_extant_tips;
    if (answ < 0.5) answ = 1.0 - answ;
    return answ;
  }
};

struct node_poly {
  std::vector<node_poly*> daughters;
  int depth;
  int dist_to_tips;
  int num_extant_tips;

  node_poly() {
    depth = 0;
    dist_to_tips = 0;
    num_extant_tips = 0;
  }

  void add_daughter(node_poly* d) {
    daughters.push_back(d);
  }

  void set_depth(int parent_depth) {
    depth = 1 + parent_depth;

    if (!daughters.empty()) {
      for (auto& i : daughters) {
        i->set_depth(depth);
      }
    }

    return;
  }

  int set_dist_to_tips() {
    dist_to_tips = 0;
    if (!daughters.empty()) {
      std::vector<double> d(daughters.size());
      for (size_t i = 0; i < daughters.size(); ++i) {
        d[i] = 1 + daughters[i]->set_dist_to_tips();
      }
      dist_to_tips = *std::max_element(d.begin(), d.end());
    }
    return dist_to_tips;
  }

  int get_acc_num_tips() {
    if (daughters.empty()) {
      num_extant_tips = 1;
    } else {
      if (daughters.size() == 1) {
        num_extant_tips = 1 + daughters.front()->get_acc_num_tips();
      } else {
        for (auto& i : daughters) {
          num_extant_tips += i->get_acc_num_tips();
        }
      }
    }

    return num_extant_tips;
  }

  bool is_double_cherry() {
    if (num_extant_tips != 4) return false;
    if (daughters.size() != 2) return false;

    int L = daughters[0]->num_extant_tips;
    int R = daughters[1]->num_extant_tips;
    if (L == 2 && R == 2) return true;

    return false;
  }

  bool is_four_prong() {
    if (num_extant_tips != 4) return false;
    if (daughters.size() != 2) return false;

    int L = daughters[0]->num_extant_tips;
    int R = daughters[1]->num_extant_tips;
    if (L == 1 && R == 3) return true;
    if (L == 3 && R == 1) return true;

    return false;
  }

  double calc_root_imbalance() {
    // this can not happe, but just to be sure:
    if (daughters.size() != 2) throw "root is not binary";

    int L = daughters[0]->num_extant_tips;

    double answ = 1.0 * L / (num_extant_tips);
    if (answ < 0.5) answ = 1.0 - answ;
    return answ;
  }
};

struct tree_base {
  virtual ~tree_base() = default;
  virtual int calc_max_width() = 0;
  virtual int calc_max_del_width() = 0;
  virtual double calc_b2() = 0;
  virtual double var_leaf_depth() = 0;
  virtual double calc_tot_int_path() = 0;
  virtual double calc_avg_vert_depth() = 0;
  virtual double calc_b1() = 0;
  virtual int max_depth() = 0;
  virtual int calc_sackin() = 0;
  virtual double calc_blum() = 0;
  virtual size_t count_cherries() = 0;
  virtual size_t count_pitchforks() = 0;
  virtual double calc_tot_coph() = 0;
  virtual size_t count_double_cherries() = 0;
  virtual size_t count_four_prong() = 0;
  virtual double root_imbalance() = 0;
};

enum used_for {depth, b1, sackin};

template<typename NODE>
class depth_tree : public tree_base {
  std::vector<NODE> tree;
  int root_no;

 public:
  explicit depth_tree(const std::vector< int >& tree_edge,
                      used_for setting = used_for::depth) {
    root_no = tree_edge[0];
    int tree_size = 0;

    for (size_t i = 2; i < tree_edge.size(); i+=2) {
      if (tree_edge[i] < root_no) root_no = tree_edge[i];
      if (tree_edge[i] > tree_size) tree_size = tree_edge[i];
    }

    tree = std::vector<NODE>(tree_size + 1);

    for (size_t i = 0; i < tree_edge.size(); i += 2) {
      auto index = static_cast<int>(tree_edge[i]);
      auto d1_index = static_cast<int>(tree_edge[i + 1]);

      tree[index].add_daughter(&tree[d1_index]);
    }
    if (setting == used_for::depth)  tree[root_no].set_depth(-1);
    if (setting == used_for::sackin) tree[root_no].get_acc_num_tips();
  }

  int calc_max_width() override {
    std::vector<int> depths(tree.size(), 0);
    for (auto i = tree.begin() + 1; i < tree.end(); ++i) {
      depths[ (*i).depth ]++;
    }
    return *std::max_element(depths.begin(), depths.end());
  }

  int calc_max_del_width() override  {
    std::vector<int> depths(tree.size(), 0);
    for (auto i = tree.begin() + 1; i < tree.end(); ++i) {
      depths[ (*i).depth ]++;
    }
    std::vector<int> dW(depths.size() - 1);
    for (size_t i = 1; i < depths.size(); ++i) {
      dW[i - 1] = depths[i] - depths[i - 1];
    }
    return(*std::max_element(dW.begin(), dW.end()));
  }

  double calc_b2() override {
    double s = 0.0;
    for (int i = 1; i < root_no; ++i) {
      // we are only interested in tip depths
      s += tree[i].depth / std::pow(2, tree[i].depth);
    }
    return s;
  }

  double var_leaf_depth() override {
    double average_depth = 0.0;
    int n = root_no - 1;
    for (int i = 1; i < root_no; ++i) {
      average_depth += tree[i].depth;
    }

    average_depth *= 1.0 / n;

    double var_depth = 0.0;
    for (int i = 1; i < root_no; ++i) {
      var_depth += (tree[i].depth - average_depth) *
        (tree[i].depth - average_depth);
    }
    var_depth *= 1.0 / n;
    return var_depth;
  }

  double calc_tot_int_path() override {
    double sum_depth = 0.0;
    for (size_t i = static_cast<size_t>(root_no); i < tree.size(); ++i) {
      sum_depth += tree[i].depth;
    }
    return sum_depth;
  }

  double calc_b1() override {
    tree[root_no].set_dist_to_tips();

    double b1 = 0.0;
    for (size_t i = static_cast<size_t>(root_no + 1); i < tree.size(); ++i) {
      b1 += 1.0 / (tree[i].dist_to_tips);
    }
    return b1;
  }

  double calc_avg_vert_depth() override {
    double sum_depth = 0.0;
    for (size_t i = 1; i < tree.size(); ++i) {
      sum_depth += tree[i].depth;
    }
    auto answ =  sum_depth * 1.0 / (tree.size() - 1);
    return answ;
  }

  int max_depth() override {
    int md = 0;
    for (const auto& i : tree) {
      if (i.depth > md) md = i.depth;
    }
    return md;
  }

  int calc_sackin() override {
    int s = 0;
    for (size_t i = root_no; i < tree.size(); ++i) {
      s += tree[i].num_extant_tips;
    }

    return s;
  }

  size_t count_pitchforks() override {
    size_t num_pitchforks = 0;
    for (size_t i = root_no; i < tree.size(); ++i) {
      if (tree[i].num_extant_tips == 3) {
        num_pitchforks++;
      }
    }
    return num_pitchforks;
  }

  size_t count_cherries() override {
    size_t num_cherries = 0;
    for (size_t i = root_no; i < tree.size(); ++i) {
      if (tree[i].num_extant_tips == 2) {
        num_cherries++;
      }
    }
    return num_cherries;
  }

   size_t count_double_cherries() override {
     size_t num_double_cherries = 0;
     for (size_t i = root_no; i < tree.size(); ++i) {
       if (tree[i].is_double_cherry()) {
         num_double_cherries++;
       }
     }
     return num_double_cherries;
   }

   size_t count_four_prong() override {
     size_t num_four_prong = 0;
     for (size_t i = root_no; i < tree.size(); ++i) {
       if (tree[i].is_four_prong()) {
         num_four_prong++;
       }
     }
     return num_four_prong;
   }

  double calc_blum() override {
    double s = 0;
    for (size_t i = root_no; i < tree.size(); ++i) {
      if (tree[i].num_extant_tips > 1) {
        s += log2(1.0 * tree[i].num_extant_tips - 1);
      }
    }
    return s;
  }

  double calc_tot_coph() override {
    double s = 0.0;
    for (size_t i = root_no + 1; i < tree.size(); ++i) {
      if (tree[i].num_extant_tips > 0) {
        s += binom_coeff(tree[i].num_extant_tips, 2);
      }
    }
    return s;
  }

  double root_imbalance() override {
    // we have previously verified the root is binary
    return tree[root_no].calc_root_imbalance();
  }
};

inline std::unique_ptr<tree_base> create_depth_tree(
    const std::vector<int>& tree_edge,
    used_for setting = used_for::depth) {
  if (is_binary(tree_edge)) {
    return std::make_unique<depth_tree<node_binary>>(tree_edge, setting);
  } else {
    return std::make_unique<depth_tree<node_poly>>(tree_edge, setting);
  }
}

}   // namespace depth

