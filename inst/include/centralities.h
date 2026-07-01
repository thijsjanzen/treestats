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
#include <numeric>
#include <algorithm>   // min element
#include <utility>     // swap
#include <stack>

#include "binom.h"  // NOLINT [build/include_subdir]

using edge   = std::vector< std::array< size_t, 2 >>;
using ltable = std::vector< std::array< double, 4>>;


struct LRSizes {

  void compute(const edge& e,
               const std::vector<double>& el,
               bool is_poly) {

    if (is_poly) {
      root_no = e[0][0];
      max_node_no = e[0][0];

      for (const auto& i : e) {
        if (i[0] < root_no) root_no = i[0];
        if (i[0] > max_node_no) max_node_no = i[0];
      }

      n_ = root_no - 1;
      N_ = max_node_no;
    } else {
      n_ = 1 + static_cast<int>(static_cast<int>(e.size()) * 0.5);  // num tips
      N_ = 2 * n_ - 1;
    }

    Tab = std::vector< std::vector< double >>(n_ - 1);
    for (auto& i : Tab) i.reserve(2); // already reserve 2 daughters

    std::array<int, 2> curRow;

    for (int ind = N_ - 2; ind >= 0; ind--) {
      curRow = { static_cast<int>(e[ind][0]) - n_ - 1,
                 static_cast<int>(e[ind][1]) - n_ - 1 };

      if (ind < 0 || ind >= static_cast<int>(el.size())) {
        throw "ind out of range el";
      }

      double W = use_branch_length ? el[ind] : 1.0;

      if (curRow[1] >= static_cast<int>(Tab.size())) {
        throw "curRow[1] out of range Tab";
      }

      double new_val = W;
      if (use_max) {
        if (curRow[1] > 0) {
          new_val += *std::max_element(Tab[curRow[1]].begin(),
                                       Tab[curRow[1]].end());
        }
      } else {
        // use sum
        if (curRow[1] > 0) {
          for (size_t i = 0; i < Tab[curRow[1]].size(); ++i) {
            new_val += Tab[curRow[1]][i];
          }
        }
      }

      int x = curRow[0];
      if (x < 0 || x >= static_cast<int>(Tab.size())) {
        throw "curRow[0] out of range Tab";
      }

      Tab[x].push_back(new_val);
    }
  }

  std::vector<double> create_q(const int offset) {
    std::vector<double> q(Tab.size(), 0.0);
    size_t cnt = 0;
    for (const auto& i : Tab) {
      q[cnt] = offset + std::accumulate(i.cbegin(), i.cend(), 0.0);
      cnt++;
    }
    return q;
  }

  LRSizes(bool use_br, bool use_m) : use_branch_length(use_br), use_max(use_m) {
    n_ = N_ = -1; // will be set later
  }

  int n_;
  int N_;
  bool use_branch_length;
  bool use_max;
  int root_no;
  int max_node_no;
  std::vector< std::vector<double> > Tab; //(n - 1, { -1, -1 });
};

double wiener(const edge& e,
              const std::vector<double>& el,
              bool normalize = false,
              bool weight = false,
              bool is_poly = false) {

  LRSizes sub_tree_sizes(weight, false);
  sub_tree_sizes.compute(e, el, is_poly);
  auto q = sub_tree_sizes.create_q(1);

  int n = static_cast<int>(q.size());
  int N = is_poly ? sub_tree_sizes.max_node_no : 2 * n + 1;


  double W = 0.0;
  for (size_t i = 0; i < e.size(); ++i) {
    int curEndpoint = static_cast<int>(e[i][1]);

    auto curQ = (curEndpoint > (n + 2)) ? q[curEndpoint - n - 2] : 1.0;

    double bl = weight ? el[i] : 1.0;

    W += curQ * (N - curQ) * bl;
  }

  if (normalize) {
    W *= 1.0 / binom_coeff_2(N);
  }

  return W;
}

double max_betweenness(const edge& e,
                       const std::vector<double>& el,
                       bool is_poly) {

  LRSizes sub_tree_sizes(false, false);

  sub_tree_sizes.compute(e, el, is_poly);

  auto q = sub_tree_sizes.create_q(0);

  auto n = q.size();

  double max_betweenness = -1.0;
  for (size_t i = 0; i < sub_tree_sizes.Tab.size(); ++i) {

    double prod = 1.0;
    for (const auto& j : sub_tree_sizes.Tab[i]) {
      prod *= j;
    }

    auto local_b = prod + q[i] * (2 * n - q[i]);
    // std::cerr << i << " " << local_b << "\n";
    // cat(i, local_b, "\n")
    if (local_b > max_betweenness) max_betweenness = local_b;
  }
  return max_betweenness;
}

double sum_weighed_heights(const edge& e,
                           const std::vector<double>& el,
                           int n = -1,
                           int N = -1) {
  if (n < 0) {
    n = 1 + static_cast<int>(static_cast<int>(e.size()) * 0.5);
    N = 2 * n - 1;
  }

  std::vector<double> Tab(N, 0.0);
  for (int ind = 0; ind < N - 1; ++ind) {
    auto curRow = e[ind];

    if (curRow[1] - 1 < 0 || curRow[1] - 1 > Tab.size()) {
      throw "curRow[1] in weighed_heights out of range";
    }
    if (curRow[0] - 1 < 0 || curRow[0] - 1 > Tab.size()) {
      throw "curRow[0] in weighed_heights out of range";
    }
    if (ind < 0 || ind >= static_cast<int>(el.size())) {
      throw "ind out of range in weighed_heights";
    }

    Tab[curRow[1] - 1] = el[ind] + Tab[curRow[0] - 1];
  }

  return std::accumulate(Tab.begin(), Tab.end(), 0.0);
}


double min_farness(const edge& local_edge,
                   const std::vector<double>& el,
                   bool weight,
                   bool is_poly) {

  LRSizes sub_tree_sizes(weight, false);
  sub_tree_sizes.compute(local_edge, el, is_poly);

  auto sizes = sub_tree_sizes.create_q(0);

  size_t n = sub_tree_sizes.n_;
  int    N = sub_tree_sizes.N_;

  std::vector<double> farness(N);

  if (n >= farness.size()) {
    throw "n >= farness.size()";
  }

  if (weight) {
    farness[n] = sum_weighed_heights(local_edge, el);
  } else {
    farness[n] = std::accumulate(sizes.begin(), sizes.end(), 0.0);
  }

  for (size_t ind = 0; ind < local_edge.size(); ++ind) {
    auto curRow = local_edge[ind];
    auto kid = curRow[1];

    if (kid > n) {
      if (kid - n - 1 < 0 || kid - n - 1 >= sizes.size()) {
        throw "kid - n - 1 outside range";
      }
    }

    double subSize = kid > n ? 1.0 + sizes[kid - n - 1] : 1.0;

    double W = weight ? el[ind] : 1.0;

    if (kid - 1 < 0 || kid - 1 >= farness.size()) {
      throw "kid outside range";
    }
    if (curRow[0] - 1 < 0 || curRow[0] - 1 >= farness.size()) {
      throw "curRow outside range";
    }

    farness[kid - 1] = farness[curRow[0] - 1] + (N - 2 * subSize) * W;
  }

  return *std::min_element(farness.begin(), farness.end());
}

double diameter(const edge& e,
                const std::vector<double>& el,
                bool weight,
                bool is_poly) {

  LRSizes depths(weight, true);
  depths.compute(e, el, is_poly);
  auto diameter = depths.create_q(0);
  return *std::max_element(diameter.begin(), diameter.end());
}

struct Edge {
  int to;
  double w;
};

void dfs(int start,
         const std::vector<std::vector<Edge>>& adj,
         std::vector<double>& dist) {
  std::fill(dist.begin(), dist.end(), -1.0);

  std::stack<int> st;
  st.push(start);
  dist[start] = 0;

  while (!st.empty()) {
    int v = st.top();
    st.pop();

    for (auto &e : adj[v]) {
      if (dist[e.to] < 0) {
        dist[e.to] = dist[v] + e.w;
        st.push(e.to);
      }
    }
  }
}

double diameter_poly(const edge& e,
                     const std::vector<double>& el,
                     bool weight = false) {
  // chatgpt magic
  int max_num_nodes = e[0][0];
  for (const auto& i : e) {
    if (i[0] > max_num_nodes) max_num_nodes = i[0];
  }

  std::vector<std::vector<Edge>> adj(max_num_nodes);

  for (int i = 0; i < e.size(); i++) {
    int a = e[i][0] - 1;
    int b = e[i][1] - 1;
    double w = weight ? el[i] : 1.0;

    adj[a].push_back({b, w});
    adj[b].push_back({a, w});
  }

  std::vector<double> dist(max_num_nodes, -1);
  dfs(0, adj, dist);

  int u = std::max_element(dist.begin(), dist.end()) - dist.begin();

  dfs(u, adj, dist);

  return *std::max_element(dist.begin(), dist.end());
}



// LTABLE associated code
class LRsizes {
public:
  explicit LRsizes(const ltable& l_in) : ltable_(l_in) {
    extant_tips = std::vector<int>(l_in.size(), 2);
    dist_to_tips = std::vector<double>(l_in.size(), 0.0);
    num_tips = get_num_tips();
  }

  std::vector<std::array<double, 2>> collect_stat_noW() {
    std::vector<std::array<double, 2>> stat;
    while (true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {  // we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);
      if (j_parent < 0) {
        throw "out of bounds";
      }

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      stat.push_back({L-1.0, R-1.0});

      if (ltable_.size() == 1) break;
    }
    return stat;
  }

  std::vector<std::array<double, 2>> collect_diameter_noW() {
    std::vector<std::array<double, 2>> stat;
    std::vector<int> depth_tips(ltable_.size(), 1);
    while (true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {    // we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = depth_tips[j];
      int R = depth_tips[j_parent];
      depth_tips[j_parent] = 1 + std::max(L, R);
      remove_from_dataset(j);

      stat.push_back({static_cast<double>(L), static_cast<double>(R)});

      if (ltable_.size() == 1) break;
    }
    return stat;
  }

  std::vector<std::array<double, 2>> collect_diameter_W() {
    std::vector<std::array<double, 2>> stat(ltable_.size() - 1);
    for (size_t i = 1; i < ltable_.size(); ++i) {
      stat[i - 1] = {ltable_[i][0], ltable_[i][0]};
    }
    return stat;
  }

  int index_of_parent(int parent) {
    int index = 0;
    bool found = false;
    for (; index < static_cast<int>(ltable_.size()); ++index) {
      if (ltable_[index][2] == parent) {
        found = true;
        break;
      }
    }
    if (!found) index = -1;
    return index;
  }

  size_t get_min_index() {
    auto min_val = std::min_element(ltable_.begin(), ltable_.end(),
                                    [&](const auto& a, const auto& b) {
                                      return a[0] < b[0];
                                    });
    return std::distance(ltable_.begin(), min_val);
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
  std::vector<double> dist_to_tips;
  std::vector<int> depth_tips;
  size_t num_tips;
};

double max_betweenness_ltable(const ltable& ltab_) {
  LRsizes left_right(ltab_);
  auto sub_tree_sizes = left_right.collect_stat_noW();

  std::vector<double> q(sub_tree_sizes.size());
  size_t cnt = 0;
  for (const auto& i : sub_tree_sizes) {
    q[cnt] = i[0] + i[1];
    cnt++;
  }
  auto n = q.size();

  double max_betweenness = -1.0;
  for (size_t i = 0; i < sub_tree_sizes.size(); ++i) {
    auto local_b = sub_tree_sizes[i][0] * sub_tree_sizes[i][1] +
      q[i] * (2 * n - q[i]);
    if (local_b > max_betweenness) max_betweenness = local_b;
  }
  return max_betweenness;
}

double diameter_ltable(const ltable& ltab_,
                       bool weight) {
  LRsizes left_right(ltab_);

  std::vector< std::array<double, 2 >> depths;

  if (weight) {
    depths = left_right.collect_diameter_W();
  } else {
    depths = left_right.collect_diameter_noW();
  }

  double diam = 0.0;
  for (const auto& i : depths) {
    auto local_depth = i[0] + i[1];
    if (local_depth > diam) diam = local_depth;
  }
  return diam;
}
