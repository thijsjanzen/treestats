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
//
//
/// other statistics


#include <vector>
#include <array>
#include <RcppArmadillo.h>

#include "util.h"        // NOLINT [build/include_subdir]
#include "mpd.h"         // NOLINT [build/include_subdir]
#include "mntd.h"        // NOLINT [build/include_subdir]

// [[Rcpp::export]]
double calc_mpd_cpp(const std::vector<int>& edge,
                     const std::vector<double>& el) {
  mpd_tree::phylo_tree focal_tree(edge, el);
  auto mpd = focal_tree.calculate_mpd();
  return mpd;
}

// [[Rcpp::export]]
double calc_J_cpp(const std::vector<int>& edge,
                  const std::vector<double>& el) {
  mpd_tree::phylo_tree focal_tree(edge, el);
  auto mpd = focal_tree.calculate_mpd();
  int n = (el.size() + 2) * 0.5;

  return mpd * 1.0 / n;
}

// [[Rcpp::export]]
double calc_var_mpd_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return calc_var_mpd_stat(edge, el);
}

// [[Rcpp::export]]
double calc_mntd_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return calc_mntd_stat(edge, el);
}

// [[Rcpp::export]]
double calc_mntd_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto ltab = convert_to_ltable(ltable_R);
  return calc_mntd_ltable(ltab);
}
