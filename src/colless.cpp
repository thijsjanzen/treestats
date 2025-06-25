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
//// COLLESS BASED or related statistics

#include <vector>
#include <array>
#include <string>
#include <RcppArmadillo.h>

#include "util.h"     // NOLINT [build/include_subdir]
#include "colless.h"  // NOLINT [build/include_subdir]
#include "ILnumber.h" // NOLINT [build/include_subdir]
#include "root_imbal.h" // NOLINT [build/include_subdir]


// [[Rcpp::export]]
double calc_colless_cpp(const std::vector<int>& parent_list,
                        std::string normalization) {
  colless_tree::colless_tree focal_tree(parent_list);
  double output = focal_tree.calc_stat(&calc_colless);

  if (normalization == "yule") {
    size_t n  = parent_list.size() / 4 + 1;
    output = correction::correct_yule(output, n);
  }
  if (normalization == "pda") {
    size_t n  = parent_list.size() / 4 + 1;
    output = correction::correct_pda(output, n);
  }
  return output;
}

// [[Rcpp::export]]
double calc_colless_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                                std::string normalization) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  double output = static_cast<double>(s.colless());

  if (normalization == "yule") {
    output = correction::correct_yule(output, s.size());
  }
  if (normalization == "pda") {
    output = correction::correct_pda(output, s.size());
  }
  return output;
}

// [[Rcpp::export]]
double calc_eWcolless_cpp(const std::vector<int>& parent_list) {
  colless_tree::colless_tree colless_tree(parent_list);
  double s = colless_tree.calc_stat(calc_ew_colless);
  return s * 1.0 / (colless_tree.size() - 1);
}

// [[Rcpp::export]]
double calc_eWcolless_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.ew_colless();
}

// [[Rcpp::export]]
size_t ILnumber_cpp(const std::vector<int>& tree_edge) {
  return calc_IL(tree_edge);
}

// [[Rcpp::export]]
size_t ILnumber_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_IL();
}

// [[Rcpp::export]]
double calc_rquartet_cpp(const std::vector<int>& tree_edge) {
  colless_tree::colless_tree focal_tree(tree_edge);
  return 3.0 * focal_tree.calc_stat(&calc_rquartet);
}

// [[Rcpp::export]]
double calc_rquartet_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return 3.0 * c.count_rquartet();
}

// [[Rcpp::export]]
double stairs_cpp(const std::vector<int>& tree_edge) {
  colless_tree::colless_tree focal_tree(tree_edge);
  double s = focal_tree.calc_stat(&calc_stairs);
  return s * 1.0 / focal_tree.size();
}

// [[Rcpp::export]]
double stairs_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_stairs();
}

// [[Rcpp::export]]
double stairs2_cpp(const std::vector<int>& tree_edge) {
  colless_tree::colless_tree focal_tree(tree_edge);
  double s = focal_tree.calc_stat(&calc_stairs2);
  return s * 1.0 / focal_tree.size();
}

// [[Rcpp::export]]
double stairs2_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_stairs2();
}

// [[Rcpp::export]]
int calc_rogers_cpp(const std::vector<int>& parent_list) {
  colless_tree::colless_tree focal_tree(parent_list);
  return focal_tree.calc_stat(&calc_rogers);
}

// [[Rcpp::export]]
double calc_rogers_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.rogers();
}

// [[Rcpp::export]]
double calc_j_one_cpp(const std::vector<int>& parent_list) {
  colless_tree::colless_tree colless_tree(parent_list);
  return colless_tree.calc_j_one(&calc_j_one);
}

// [[Rcpp::export]]
double calc_j_one_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.collect_j_one();
}

// [[Rcpp::export]]
double calc_Ibased_cpp(const std::vector<int>& parent_list) {
  colless_tree::colless_tree colless_tree(parent_list);
  return colless_tree.collect_I(&calc_I);
}

// [[Rcpp::export]]
double calc_Ibased_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.collect_I();
}

// [[Rcpp::export]]
size_t pitchforks_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  // non-ltable version uses sackin
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_pitchforks();
}

// [[Rcpp::export]]
double calc_colless_corr_cpp(const std::vector<int>& parent_list,
                             std::string normalization) {
  colless_tree::colless_tree focal_tree(parent_list);
  double output = focal_tree.calc_stat(&calc_colless);

  size_t n = focal_tree.size() + 1;

  output *= 2.0 / ((n - 1) * (n - 2));

  if (normalization == "yule") {
    auto expected_val_yule = 2.0 * n / ((n - 1) * (n - 2));

    double sum_thing = 0.0;
    if (n % 2 == 0) {  // even number of tips
      for (size_t j = 2; j < n / 2; ++j) sum_thing += 1.0 / j;
    } else {
      sum_thing = 1.0 / n;
      for (size_t j = 2; j < (n - 1) / 2; ++j) sum_thing += 1.0 / j;
    }

    expected_val_yule *= sum_thing;

    output *= 1.0 / expected_val_yule;
  }

  return output;
}

// [[Rcpp::export]]
double calc_colless_corr_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                               std::string normalization) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  double output = static_cast<double>(s.colless());

  size_t n = s.size();
  output *= 2.0 / ((n - 1) * (n - 2));

  if (normalization == "yule") {
    auto expected_val_yule = 2.0 * n / ((n - 1) * (n - 2));

    double sum_thing = 0.0;
    if (n % 2 == 0) {  // even number of tips
      for (size_t j = 2; j < n / 2; ++j) sum_thing += 1.0 / j;
    } else {
      sum_thing = 1.0 / n;
      for (size_t j = 2; j < (n - 1) / 2; ++j) sum_thing += 1.0 / j;
    }
    expected_val_yule *= sum_thing;

    output *= 1.0 / expected_val_yule;
  }
  return output;
}

// [[Rcpp::export]]
double calc_colless_quad_cpp(const std::vector<int>& parent_list,
                             std::string normalization) {
  colless_tree::colless_tree focal_tree(parent_list);
  double output = focal_tree.calc_stat(&calc_colless_quad);
  if (normalization == "yule") {
    size_t n = focal_tree.size() + 1;
    auto expected_yule = n * (n + 1);
    double sum_thing = 0.0;
    for (size_t i = 1; i <= n; ++i) sum_thing += 1.0 / i;

    expected_yule -= 2 * n * sum_thing;

    output *= 1.0 / expected_yule;
  }

  return output;
}

// [[Rcpp::export]]
double calc_colless_quad_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                                    std::string normalization) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  double output = static_cast<double>(s.colless_quad());

  if (normalization == "yule") {
    size_t n = s.size();
    auto expected_yule = n * (n + 1);
    double sum_thing = 0.0;
    for (size_t i = 1; i <= n; ++i) sum_thing += 1.0 / i;

    expected_yule -= 2 * n * sum_thing;

    output *= 1.0 / expected_yule;
  }
  return output;
}

// [[Rcpp::export]]
double calc_root_imbalance_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  auto imbal = calc_root_imbal(l_in_cpp);
  return imbal;
}

// [[Rcpp::export]]
double calc_root_imbalance_cpp(const std::vector<int>& parent_list) {
  colless_tree::colless_tree focal_tree(parent_list);
  return focal_tree.calc_root_imbal();
}

// [[Rcpp::export]]
double calc_double_cherries_cpp(const std::vector<int>& parent_list) {
  colless_tree::colless_tree focal_tree(parent_list);
  return focal_tree.calc_double_cherries();
}

// [[Rcpp::export]]
double calc_double_cherries_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_double_cherries();
}

// [[Rcpp::export]]
double calc_four_prong_cpp(const std::vector<int>& parent_list) {
  colless_tree::colless_tree focal_tree(parent_list);
  return focal_tree.calc_four_prong();
}

// [[Rcpp::export]]
double calc_four_prong_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_four_prong();
}

