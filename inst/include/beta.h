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

#include <random>
#include <cmath>
#include <array>
#include <unordered_map>
#include <nloptrAPI.h>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <RcppArmadillo.h>


using ltable = std::vector< std::array<double, 4>>;

class betastat {
 public:
  explicit betastat(const std::vector< std::array< int, 2 >> e) : edge(e) {
    tiplist = std::vector<int>(edge.size() + 2, -1);
    update_lr_matrix();
  }

  explicit betastat(const ltable& lt_in) : lt_(lt_in) {
    for (auto i : lt_) {
      brts_.push_back(i[0]);
    }
    std::sort(brts_.begin(), brts_.end());
    brts_.erase(std::unique(brts_.begin(), brts_.end()),
                 brts_.end());
    update_lr_matrix_ltable();
  }

  double calc_likelihood(double beta) const {
    std::vector< double > sn = get_sn(beta);
    std::vector< double > ll(lr_.size());
    assert(ll.size() == n_.size());
    for (size_t i = 0; i < n_.size(); ++i) {
      auto index = n_[i];
      assert(index < sn.size() );
      ll[i] = calc_log_prob(i, sn[ index ], beta);
    }
    double sumll = std::accumulate(ll.begin(), ll.end(), 0.f);
    return sumll;
  }

 private:
  std::vector< std::array<int, 2>> lr_;
  std::vector< std::array<int, 2>> edge;
  int max_n_;
  std::vector< int > n_;

  std::vector< int > tiplist;

  const ltable lt_;
  std::vector<double> brts_;

  int get_num_tips(const int& label,
                   const int& root_label) {
    if (label >= static_cast<int>(tiplist.size()) ||
        label < 0) {
      throw std::out_of_range("label outside tiplist.size()");
    }

    if (label < root_label) {
      tiplist[label] = 1;
      return 1;
    }

    if (tiplist[label] > 0) {   // tiplist is populated with -1
      return(tiplist[label]);
    }

    std::vector< int > matches(2);
    auto match1 = std::lower_bound(edge.begin(), edge.end(), label,
                                   [&](const auto& a, int val){
      return a[0] < val;
    });

    if (match1 != edge.end()) {
      if ((*match1)[0] == label) {
        matches[0] = (*match1)[1];
        match1++;
        if ((*match1)[0] == label) {
          matches[1] = (*match1)[1];
        } else {
          matches.pop_back();
        }
      }
    } else {
      throw "can't find matches";
    }

    int s = 0;
    for (const auto& i : matches) {
      s += get_num_tips(i, root_label);
    }
    tiplist[label] = s;
    return s;
  }


  void update_lr_matrix() {
    auto root_label = edge[0][0];

    std::sort(edge.begin(), edge.end(), [&](const auto& a, const auto& b) {
      return a[0] < b[0];
    });

    for (size_t i = 0; i < edge.size(); ++i) {
      if (i + 1 < edge.size()) {
        auto j = i + 1;
        std::array< int, 2> lr;
        lr[0] = get_num_tips(edge[i][1], root_label);
        lr[1] = get_num_tips(edge[j][1], root_label);

        if (lr[0] > lr[1]) {
          std::swap(lr[0], lr[1]);
        }
        size_t total_num_lin = lr[0] + lr[1];
        n_.push_back(total_num_lin);
        lr_.push_back(lr);
        i++;
      }
    }

    max_n_ = *std::max_element(n_.cbegin(), n_.cend());
  }


  double calc_i_n_b(int i, int n, double b) const {
    double nom   = std::tgamma(1.f*(i + 1 + b)) *
                   std::tgamma(1.f*(n - i + 1 + b));
    double denom = std::tgamma(1.f*(i + 1)) *
                   std::tgamma(1.f*(n - i + 1));
    return(nom / denom);
  }

  double calc_i_n_b_l(int i, int n, double b) const {
    return lgamma(i + 1 + b) + lgamma(n - i + 1 + b) -
      lgamma(i + 1) - lgamma(n - i + 1);
  }

  std::vector<double> get_sn(double b) const {
    std::vector<double> sn(max_n_ + 1, 0.0);
    std::vector<double> xn(max_n_ + 1, 0.0);

    if (sn.size() < 4) {
      throw std::out_of_range("get_n too small tree");
    }

    xn[2] = 1.0;
    xn[3] = 0.5;

    sn[2] = expf(calc_i_n_b_l(1, 2, b));
    sn[3] = expf(calc_i_n_b_l(1, 3, b)) + expf(calc_i_n_b_l(2, 3, b));

    if (max_n_ >= 3) {
      for (int n = 3; n < max_n_; ++n) {
        auto term1 = n + 2 + 2 * b;
        auto term2 = 2 * (n + b) * xn[n];
        xn[n + 1] = ((n + b) * (n + 1) * xn[n]) / (n * term1 + term2);
        sn[n + 1] =  (1.0 / (n + 1) ) * (term1 + term2 / n) * sn[n];
      }
    }

    return sn;
  }

  double calc_log_prob(int index, double sn, double beta) const {
    double l = lr_[index][0];
    double r = lr_[index][1];
    return lgamma(beta + l + 1) + lgamma(beta + r + 1) -
           lgamma(l + 1) - lgamma(r + 1) - log(sn);
  }

  int find_species_in_ltable(int sp) {
    for (int i = 0; i < static_cast<int>(lt_.size()); ++i) {
      if (lt_[i][2] == sp) {
        return i;
      }
    }
    throw "can't find species in ltable\n";
    return -1;
  }

  std::vector<double> find_daughters(int sp,
                                     double bt) {
    std::vector< double > output;
    for (auto i : lt_) {
      if (i[0] < bt && i[1] == sp) {
        output.push_back(i[2]);
      }
    }
    return(output);
  }

  int get_total_num_lin(int sp,
                           double bt) {
    int index = find_species_in_ltable(sp);
    int total_tips = 0;
    if (index >= 0) {
      if (lt_[index][3] == -1) {
        total_tips = 1;
      }
    }

    // now, we have to find daughters branching off
    std::vector< double > daughters = find_daughters(sp, bt);
    if (!daughters.empty())  {
      for (auto d : daughters) {
        total_tips += get_total_num_lin(static_cast<int>(d), bt);
      }
    }
    return(total_tips);
  }

  std::vector< size_t > get_indices(double bt) {
    std::vector< size_t > indices;
    for (size_t i = 0; i < lt_.size(); ++i) {
      if (lt_[i][0] == bt) {
        indices.push_back(i);
      }
    }
    return indices;
  }

  /// ltable member functions
  void update_lr_matrix_ltable() {
    max_n_ = 0;
    for (auto br : brts_) {
      std::array<int, 2> lr = {0, 0};
      std::vector< size_t > indices = get_indices(br);
      if (indices.size() == 2) {
        lr[0] = get_total_num_lin(lt_[indices[0]][2], br);
        lr[1] = get_total_num_lin(lt_[indices[1]][2], br);
      }
      if (indices.size() == 1) {
        lr[0] = get_total_num_lin(lt_[indices[0]][2], br);
        lr[1] = get_total_num_lin(lt_[indices[0]][1], br);
      }

      if (lr[0] > lr[1]) {
        std::swap(lr[0], lr[1]);
      }
      int total_num_lin = lr[0] + lr[1];
      if (total_num_lin > max_n_) max_n_ = total_num_lin;
      n_.push_back(total_num_lin);
      lr_.push_back(lr);
    }

    return;
  }
};

struct nlopt_f_data {
  explicit nlopt_f_data(const betastat& b_in) : b(b_in) {
  }
  const betastat b;
};

double objective(unsigned int n, const double* x, double*, void* func_data) {
  auto psd = reinterpret_cast<nlopt_f_data*>(func_data);
  return(-psd->b.calc_likelihood(x[0]));
}

template< class T>
double calc_beta(const T& edge,
                 double lower_lim,
                 double upper_lim,
                 std::string algorithm,
                 double abs_tol,
                 double rel_tol) {
  betastat beta_calc(edge);
  // now we do optimization

  nlopt_f_data optim_data(beta_calc);

  nlopt_opt opt;
  bool algo_set = false;
  double init_val = -1.9;

  if (algorithm == "subplex") {
    opt = nlopt_create(NLOPT_LN_SBPLX, static_cast<unsigned int>(1));
    algo_set = true;
  }
  if (algorithm == "simplex") {
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, static_cast<unsigned int>(1));
    algo_set = true;
  }
  if (algorithm == "COBYLA") {
    opt = nlopt_create(NLOPT_LN_COBYLA, static_cast<unsigned int>(1));
    algo_set = true;
    init_val = 0.01;
  }

  if (!algo_set) {
     throw "no algorithm chosen";
  }

  double llim[1] = {static_cast<double>(lower_lim)};
  double ulim[1] = {static_cast<double>(upper_lim)};

  nlopt_set_lower_bounds(opt, llim);
  nlopt_set_upper_bounds(opt, ulim);

  nlopt_set_min_objective(opt, objective, &optim_data);

  nlopt_set_xtol_rel(opt, rel_tol);
  nlopt_set_ftol_abs(opt, abs_tol);

  std::vector<double> x = {init_val};
  double minf;

  auto nloptresult = nlopt_optimize(opt, &(x[0]), &minf);

  if (nloptresult < 0) {
    Rcpp::Rcout << "failure to optimize!\n";
  }

  nlopt_destroy(opt);

  return x[0];  // beta
}
