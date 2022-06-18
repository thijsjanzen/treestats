#ifndef CENTRALITIES_H
#define CENTRALITIES_H

#include <vector>
#include <array>
#include <numeric>

#include "binom.h"

using edge = std::vector< std::array< size_t, 2 >>;

// first sum:
std::vector< std::array< double, 2 >> computeLRSizes(const edge& e,
                                                     const std::vector<double>& el,
                                                     bool use_branch_length = false,
                                                     bool use_max = false) {

  int n = 1 + static_cast<int>(static_cast<int>(e.size()) * 0.5); // num tips
  int N = 2 * n - 1;

  std::vector< std::array< double, 2 >> Tab(n - 1, { -1, -1 });

  std::array<int, 2> curRow;

  for (int ind = N - 2; ind >= 0; ind--) {
    curRow = { static_cast<int>(e[ind][0]) - n - 1,
               static_cast<int>(e[ind][1]) - n - 1 };


    if (ind < 0 || ind >= el.size()) {
      throw "ind out of range el";
    }

    double W = use_branch_length ? el[ind] : 1.0;


    if (curRow[1] >= static_cast<int>(Tab.size())) {
      throw "curRow[1] out of range Tab";
    }

    double new_val;
    if (use_max) {
      new_val = curRow[1] > 0 ? W + std::max(Tab[curRow[1]][0], Tab[curRow[1]][1]) : W;
    } else {
      // use sum
      new_val = curRow[1] > 0 ? W + Tab[curRow[1]][0] + Tab[curRow[1]][1] : W;
    }

    int x = curRow[0];
    if (curRow[0] < 0 || curRow[0] >= Tab.size()) {
      throw "curRow[0] out of range Tab";
    }
    int y = Tab[curRow[0]][0] < 0 ? 0 : 1;

    Tab[x][y] = new_val;
  }

  return Tab;
}


double wiener(const edge& e,
              const std::vector<double>& el,
              bool normalize = false,
              bool weight = false) {

  auto sub_tree_sizes = computeLRSizes(e, el);
  std::vector<double> q(sub_tree_sizes.size());
  size_t cnt = 0;
  for (const auto& i : sub_tree_sizes) {
    q[cnt] = i[0] + i[1] + 1;
    cnt++;
  }
  int n = static_cast<int>(q.size());
  int N = 2 * n + 1;

  double W;
  if (weight) {
    W = 0.0;
    for (size_t i = 0; i < e.size(); ++i) {
      int curEndpoint = static_cast<int>(e[i][1]);

      auto curQ = (curEndpoint > (n + 2)) ? q[curEndpoint - n - 2] : 1.0;

      W += curQ * (N - curQ) * el[i];
    }
  }
  else {
    W = (N - 1) * (n + 1);
    for (const auto& i : q) {
      W += i * (N - i);
    }
  }

  if (normalize) {
    W *= 1.0 / binom_coeff_2(N);
  }

  return W;
}


double max_betweenness(const edge& e,
                       const std::vector<double>& el) {

  auto sub_tree_sizes = computeLRSizes(e, el);
  std::vector<double> q(sub_tree_sizes.size());
  size_t cnt = 0;
  for (const auto& i : sub_tree_sizes) {
    q[cnt] = i[0] + i[1];
    cnt++;
  }
  auto n = q.size();

  double max_betweenness = -1.0;
  for (size_t i = 0; i < sub_tree_sizes.size(); ++i) {
    auto local_b = sub_tree_sizes[i][0] * sub_tree_sizes[i][1] + q[i] * (2 * n - q[i]);
    if (local_b > max_betweenness) max_betweenness = local_b;
  }
  return max_betweenness;
}

double sum_weighed_heights(const edge& e,
                           const std::vector<double>& el) {
  int n = 1 + static_cast<int>(static_cast<int>(e.size()) * 0.5);
  int N = 2 * n - 1;
  std::vector<double> Tab(N, 0.0);
  for (int ind = 0; ind < N - 1; ++ind) {
    auto curRow = e[ind];

    if (curRow[1] - 1 < 0 || curRow[1] - 1 > Tab.size()) {
      throw "curRow[1] in weighed_heights out of range";
    }
    if (curRow[0] - 1 < 0 || curRow[0] - 1 > Tab.size()) {
      throw "curRow[0] in weighed_heights out of range";
    }
    if (ind < 0 || ind >= el.size()) {
      throw "ind out of range in weighed_heights";
    }

    Tab[curRow[1] - 1] = el[ind] + Tab[curRow[0] - 1];
  }

  return std::accumulate(Tab.begin(), Tab.end(), 0.0);
}

double min_farness(const edge& local_edge,
                   const std::vector<double>& el,
                   bool weight = false) {
  auto sub_tree_sizes = computeLRSizes(local_edge, el);

 // std::cerr << "LRSizes done\n"; std::this_thread::sleep_for(std::chrono::milliseconds(300));

  std::vector<double> sizes(sub_tree_sizes.size());
  size_t cnt = 0;
  for (const auto& i : sub_tree_sizes) {
    sizes[cnt] = i[0] + i[1];
    cnt++;
  }

 // std::cerr << "sizes collected\n"; std::this_thread::sleep_for(std::chrono::milliseconds(300));


  int n = 1 + static_cast<int>(local_edge.size()) * 0.5;;
  int N = 2 * n - 1;

  std::vector<double> farness(N);

  if (n >= farness.size()) {
    throw "n >= farness.size()";
  }

  if (weight) {
    farness[n] = sum_weighed_heights(local_edge, el);
  }
  else {
    farness[n] = std::accumulate(sizes.begin(), sizes.end(), 0.0);
  }

//  std::cerr << "starting ind loop\n"; std::this_thread::sleep_for(std::chrono::milliseconds(300));
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

//  std::cerr << "starting min_element\n"; std::this_thread::sleep_for(std::chrono::milliseconds(300));

  return *std::min_element(farness.begin(), farness.end());
}

double max_closeness(const edge& e,
                     const std::vector<double>& el,
                     bool weight = false) {
  return 1.0 / min_farness(e, el, weight);
}

double diameter(const edge& e,
                const std::vector<double>& el,
                bool weight = false) {
  auto depths = computeLRSizes(e, el, weight, true);
  double diam = 0.0;
  for (const auto& i : depths) {
    auto local_depth = i[0] + i[1];
    if (local_depth > diam) diam = local_depth;
  }
  return diam;
}




#endif /* CENTRALITIES_H */
