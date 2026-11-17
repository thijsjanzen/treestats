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
    output = correction_colless::correct_yule(output, n);
  }
  if (normalization == "pda") {
    size_t n  = parent_list.size() / 4 + 1;
    output = correction_colless::correct_pda(output, n);
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
    output = correction_colless::correct_yule(output, s.size());
  }
  if (normalization == "pda") {
    output = correction_colless::correct_pda(output, s.size());
  }
  return output;
}

double normalize_ew_colless(double stat, size_t n) {
  bool is_even = !(n % 2);
  double correction = 0.0;
  if (is_even) {
    auto numerator = 4 * n * (n - 1) * (harmonic(n) - harmonic(n / 2)) - n;
    auto denumerator = 6 * (n - 1) * (n - 2);

    correction = numerator / denumerator;
  } else {
    auto numerator = 4 * n * (n - 2) * (harmonic(n - 1) -
                     harmonic((n - 1) / 2)) - n + 3;
    auto denumerator = 6 * (n - 2) * (n - 2);
    correction = numerator / denumerator;
  }
  return stat / correction;
}

// [[Rcpp::export]]
double calc_eWcolless_cpp(const std::vector<int>& parent_list,
                          std::string normalize) {
  colless_tree::colless_tree colless_tree(parent_list);
  double s = colless_tree.calc_stat(calc_ew_colless);
  auto answ = s * 1.0 / (colless_tree.size() - 1);
  if (normalize == "yule") {
    answ = normalize_ew_colless(answ, colless_tree.size());
  }
  return answ;
}

// [[Rcpp::export]]
double calc_eWcolless_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                                 std::string normalize) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  auto answ = s.ew_colless();
  if (normalize == "yule") {
    answ = normalize_ew_colless(answ, l_in_cpp.size());
  }
  return answ;
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

double correct_stairs1(double stat, size_t n) {
  // correction under the expectation of the Yule model
  bool is_even = n % 2 == 0;

  double correction = is_even ?
    2 * n * (harmonic((n / 2) - 1) - harmonic(n - 1) + 1) :       // even
    2 * n * (harmonic((n - 1) / 2) - harmonic(n - 1) + 1) - 2;    // odd

  correction *= 1.0 / (n - 1);

  stat *= 1.0 / correction;
  return stat;
}

double correct_stairs2(double stat, size_t n) {
  double correction;
  if (n % 2 == 0) {     // even
    correction = 2 * n * harmonic2(n - 1) - n * harmonic2(n / 2)  -
                 2 * (harmonic(n - 1) - harmonic(n / 2)) - (n - 1);
  } else {              // odd
    correction = 2 * n * harmonic2(n - 1) - n * harmonic2((n - 1) / 2)  -
                 2 * (harmonic(n - 1) - harmonic((n - 1) / 2)) - (n - 1);
  }
  correction *= 1.0 / (n - 1);
  stat *= 1.0 / correction;
  return stat;
}

// [[Rcpp::export]]
double stairs_cpp(const std::vector<int>& tree_edge,
                  std::string normalization) {
  colless_tree::colless_tree focal_tree(tree_edge);
  double s = focal_tree.calc_stat(&calc_stairs);

  auto answ = s * 1.0 / focal_tree.size();
  if (normalization == "yule") {
    size_t n  = tree_edge.size() / 4 + 1;
    answ = correct_stairs1(answ, n);
  }
  return answ;
}

// [[Rcpp::export]]
double stairs_ltable_cpp(const Rcpp::NumericMatrix& ltable_R,
                         std::string normalization) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  auto answ = c.count_stairs();
  if (normalization == "yule") {
    answ = correct_stairs1(answ, local_ltab.size());
  }
  return answ;
}

// [[Rcpp::export]]
double stairs2_cpp(const std::vector<int>& tree_edge,
                   std::string normalization) {
  colless_tree::colless_tree focal_tree(tree_edge);
  double s = focal_tree.calc_stat(&calc_stairs2);
  auto answ = s * 1.0 / focal_tree.size();
  if (normalization == "yule") {
    answ = correct_stairs2(answ, focal_tree.size());
  }
  return answ;
}

// [[Rcpp::export]]
double stairs2_ltable_cpp(const Rcpp::NumericMatrix& ltable_R,
                          std::string normalization) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  auto answ = c.count_stairs2();
  if (normalization == "yule") {
    answ = correct_stairs2(answ, local_ltab.size());
  }
  return answ;
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

double colless_corr_yule(size_t n) {
  double numerator = 0.0;
  if (n % 2 == 0) {  // even number of tips
    numerator = 2 * n * (harmonic(n / 2) - 1);
  } else {
    numerator = 2 * n * (harmonic((n - 1) / 2) - 1) + 2;
  }
  auto answer = numerator * 1.0 / ((n - 2) * (n - 1));
  return answer;
}

// [[Rcpp::export]]
double calc_colless_corr_cpp(const std::vector<int>& parent_list,
                             std::string normalization) {
  colless_tree::colless_tree focal_tree(parent_list);
  double output = focal_tree.calc_stat(&calc_colless);

  size_t n = focal_tree.size() + 1;

  output *= 2.0 / ((n - 1) * (n - 2));

  if (normalization == "yule") {
    auto expected_val_yule = colless_corr_yule(n);

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
    double expected_val_yule = colless_corr_yule(n);

    output *= 1.0 / expected_val_yule;
  }
  return output;
}

double colless_corr_quad(size_t n) {
  auto expected_yule = n * (n + 1);
  double sum_thing = 0.0;
  for (size_t i = 1; i <= n; ++i) sum_thing += 1.0 / i;

  expected_yule -= 2 * n * sum_thing;
  return expected_yule;
}

// [[Rcpp::export]]
double calc_colless_quad_cpp(const std::vector<int>& parent_list,
                             std::string normalization) {
  colless_tree::colless_tree focal_tree(parent_list);
  double output = focal_tree.calc_stat(&calc_colless_quad);
  if (normalization == "yule") {
    size_t n = focal_tree.size() + 1;
    auto expected_yule = colless_corr_quad(n);

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
    auto expected_yule = colless_corr_quad(n);
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
double calc_double_cherries_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_double_cherries();
}

// [[Rcpp::export]]
double calc_four_prong_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_four_prong();
}
