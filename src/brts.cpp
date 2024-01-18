// Copyright 2022 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//
// /// BRTS based statistics

#include <vector>
#include <array>
#include <Rcpp.h>

#include "branching_times.h" // NOLINT [build/include_subdir]
#include "phylo2L.h"
#include "pigot_rho.h"       // NOLINT [build/include_subdir]
#include "gamma.h"           // NOLINT [build/include_subdir]
#include "nltt.h"            // NOLINT [build/include_subdir]
#include "crown_age.h"       // NOLINT [build/include_subdir]

#include "util.h"            // NOLINT [build/include_subdir]


// [[Rcpp::export]]
std::vector< double > branching_times_cpp(const Rcpp::List& phy) {
  return branching_times_phy(phy);
}

// [[Rcpp::export]]
std::vector<double>
  branching_times_ltable_cpp(const Rcpp::NumericMatrix& mat_in) {
  std::vector<double> out(mat_in.nrow() - 1);
  for (int i = 1; i < mat_in.nrow(); ++i) {
    out[i - 1] = mat_in(i, 0);
  }
  return out;
}

// [[Rcpp::export]]
double calc_rho_complete_cpp(const Rcpp::List& phy) {
  Rcpp::NumericMatrix edge = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  std::vector<double> el(edge_length.begin(), edge_length.end());
  edge_table edges(edge.nrow());
  for (int i = 0; i < edge.nrow(); i++) {
    std::array<size_t, 2> to_add = {static_cast<size_t>(edge(i, 0)),
                                    static_cast<size_t>(edge(i, 1))};
    edges[i] = to_add;
  }

  double crown_age = calc_crown_age(edges, el);
  phylo phylo_tree(edges, el);
  rho pigot_rho(phylo_tree, crown_age);
  return pigot_rho.calc_pigot_rho();
}


// [[Rcpp::export]]
double calc_rho_cpp(const Rcpp::List& phy) {
  size_t num_nodes = static_cast<size_t>(phy["Nnode"]);

  if (num_nodes < 200) {
    return calc_rho_complete_cpp(phy);
  }

  auto brts = branching_times_cpp(phy);
  return calc_rho(brts);
}

// [[Rcpp::export]]
double calc_rho_ltable_cpp(const Rcpp::NumericMatrix& ltab) {
  auto brts = branching_times_ltable_cpp(ltab);
  return calc_rho(brts);
}

// [[Rcpp::export]]
double calc_phylodiv_cpp(const Rcpp::List& phy,
                         double t,
                         double extinct_acc) {
  try {
    Rcpp::NumericMatrix edge = phy["edge"];
    Rcpp::NumericVector edge_length = phy["edge.length"];

    std::vector<double> el(edge_length.begin(), edge_length.end());
    edge_table edges(edge.nrow());
    for (int i = 0; i < edge.nrow(); i++) {
      std::array<size_t, 2> to_add = {static_cast<size_t>(edge(i, 0)),
                                   static_cast<size_t>(edge(i, 1))};
      edges[i] = to_add;
    }

    double crown_age = calc_crown_age(edges, el);  // ignore root edge
    phylo phylo_tree(edges, el);

    // function below calculates [0, T], max_t is in [T, 0]
    t = crown_age - t;

    return calculate_phylogenetic_diversity(phylo_tree,
                                            t,
                                            crown_age,
                                            extinct_acc);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}


// [[Rcpp::export]]
double calc_phylodiv_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  return calculate_phy_div_ltable(local_ltab);
}

// [[Rcpp::export]]
double calc_mean_branch_length_cpp(const Rcpp::NumericVector& edge_length) {
  auto sum_bl = std::accumulate(edge_length.begin(), edge_length.end(), 0.0);
  return sum_bl * 1.0 / (edge_length.size());
}

// [[Rcpp::export]]
double calc_mean_branch_length_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  double sum_bl = calculate_phy_div_ltable(local_ltab);
  double n = ltable_R.nrow();
  double num_branches = n + n - 2;
  return sum_bl * 1.0 / num_branches;
}

// [[Rcpp::export]]
double calc_gamma_cpp(const Rcpp::List& phy) {
  std::vector<double> brts = branching_times_cpp(phy);
  return calc_gamma(brts);
}

// [[Rcpp::export]]
double calc_gamma_cpp2(const std::vector<int>& edge,
                       const std::vector<double>& el) {
  return calc_gamma2(edge, el);
}


// [[Rcpp::export]]
double calc_gamma_ltable_cpp(const Rcpp::NumericMatrix& ltab_in) {
  std::vector<double> brts = branching_times_ltable_cpp(ltab_in);
  return calc_gamma(brts);
}


// [[Rcpp::export]]
double calc_nltt_cpp(const Rcpp::List& phy1,
                     const Rcpp::List& phy2) {
  std::vector<double> brts_one = branching_times_cpp(phy1);
  std::vector<double> brts_two = branching_times_cpp(phy2);
  std::sort(brts_one.begin(), brts_one.end(), std::greater<double>());
  std::sort(brts_two.begin(), brts_two.end(), std::greater<double>());
  for (auto& i : brts_one) {
    i *= -1;
  }
  for (auto& i : brts_two) {
    i *= -1;
  }
  brts_one.push_back(0.0);
  brts_two.push_back(0.0);
  auto nltt = calc_nltt(brts_one, brts_two);
  return nltt;
}


// [[Rcpp::export]]
double calc_nltt_ltable_cpp(const Rcpp::NumericMatrix& ltab1,
                            const Rcpp::NumericMatrix& ltab2) {
  auto brts_one = branching_times_ltable_cpp(ltab1);
  auto brts_two = branching_times_ltable_cpp(ltab2);
  std::sort(brts_one.begin(), brts_one.end(), std::greater<double>());
  std::sort(brts_two.begin(), brts_two.end(), std::greater<double>());
  for (auto& i : brts_one) {
    i *= -1;
  }
  for (auto& i : brts_two) {
    i *= -1;
  }
  brts_one.push_back(0.0);
  brts_two.push_back(0.0);
  auto nltt = calc_nltt(brts_one, brts_two);
  return nltt;
}

// [[Rcpp::export]]
double calc_crown_age_cpp(const Rcpp::List& phy) {
  Rcpp::NumericMatrix edge = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  std::vector<double> el(edge_length.begin(), edge_length.end());
  edge_table edges(edge.nrow());
  for (int i = 0; i < edge.nrow(); i++) {
    std::array<size_t, 2> to_add = {static_cast<size_t>(edge(i, 0)),
                                    static_cast<size_t>(edge(i, 1))};
    edges[i] = to_add;
  }
  return calc_crown_age(edges, el);
}

