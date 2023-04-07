// Copyright 2022 - 2023 Thijs Janzen
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
#include <Rcpp.h>

#include "util.h"     // NOLINT [build/include_subdir]
#include "colless.h"  // NOLINT [build/include_subdir]
#include "ILnumber.h" // NOLINT [build/include_subdir]


// [[Rcpp::export]]
double calc_colless_cpp(const std::vector<long>& parent_list,
                        std::string normalization) {
  colless_tree::phylo_tree colless_tree(parent_list);
  double output = static_cast<double>(colless_tree.calc_colless());

  if (normalization == "yule") {
    size_t n  = parent_list.size() / 4 + 1;
    output = colless_tree.correct_yule(output, n);
  }
  if (normalization == "pda") {
    size_t n  = parent_list.size() / 4 + 1;
    output = colless_tree.correct_pda(output, n);
  }
  return output;
}

// [[Rcpp::export]]
double calc_colless_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                                std::string normalization) {

  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  double output = static_cast<double>(s.calc_colless());

  if (normalization == "yule") {
    output = s.correct_yule(output);
  }
  if (normalization == "pda") {
    output = s.correct_pda(output);
  }
  return output;
}

// [[Rcpp::export]]
double calc_eWcolless_cpp(const std::vector<long>& parent_list) {
  colless_tree::phylo_tree colless_tree(parent_list);
  return colless_tree.calc_eWcolless();
}


// [[Rcpp::export]]
double calc_eWcolless_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_ew_colless();
}

// [[Rcpp::export]]
size_t ILnumber_cpp(const std::vector<long>& tree_edge) {
  return calc_IL(tree_edge);
}

// [[Rcpp::export]]
size_t ILnumber_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_IL();
}


// [[Rcpp::export]]
double calc_rquartet_cpp(const std::vector<long>& tree_edge,
                         std::string normalization) {
  colless_tree::phylo_tree tree(tree_edge);
  auto output = tree.calc_rquartet();

  if (normalization == "yule") {
    size_t n  = tree_edge.size() / 4 + 1;
    output = tree.correct_rquartet_yule(output, n);
  }
  if (normalization == "pda") {
    size_t n  = tree_edge.size() / 4 + 1;
    output = tree.correct_rquartet_pda(output, n);
  }
  return output;
}

// [[Rcpp::export]]
double calc_rquartet_ltable_cpp(const Rcpp::NumericMatrix& ltable_R,
                                std::string normalization) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  auto output = c.count_rquartet();

  if (normalization == "yule") {
    output = c.correct_rquartet_yule(output);
  }
  if (normalization == "pda") {
    output = c.correct_rquartet_pda(output);
  }
  return output;
}

// [[Rcpp::export]]
double stairs_cpp(const std::vector<long>& tree_edge) {
  colless_tree::phylo_tree tree(tree_edge);
  return tree.calc_stairs();
}

// [[Rcpp::export]]
double stairs_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_stairs();
}

// [[Rcpp::export]]
double stairs2_cpp(const std::vector<long>& tree_edge) {
  colless_tree::phylo_tree tree(tree_edge);
  return tree.calc_stairs2();
}

// [[Rcpp::export]]
double stairs2_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_stairs2();
}


// [[Rcpp::export]]
int calc_rogers_cpp(const std::vector<long>& parent_list) {
  colless_tree::phylo_tree colless_tree(parent_list);
  return colless_tree.calc_rogers();
}

// [[Rcpp::export]]
double calc_rogers_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_rogers();
}

// [[Rcpp::export]]
double calc_j_one_cpp(const std::vector<long>& parent_list) {
  colless_tree::phylo_tree colless_tree(parent_list);
  return colless_tree.calc_j_one();
}

// [[Rcpp::export]]
double calc_j_one_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.collect_j_one();
}


// [[Rcpp::export]]
double calc_Ibased_cpp(const std::vector<long>& parent_list) {
  colless_tree::phylo_tree colless_tree(parent_list);
  auto I_vec = colless_tree.collect_I();
  auto sum = std::accumulate(I_vec.begin(), I_vec.end(), 0.0);
  return sum * 1.0 / I_vec.size();
}

// [[Rcpp::export]]
double calc_Ibased_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  auto I_vec = s.collect_I();
  auto sum = std::accumulate(I_vec.begin(), I_vec.end(), 0.0);
  return sum * 1.0 / I_vec.size();
}

// [[Rcpp::export]]
size_t pitchforks_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  // non-ltable version uses sackin
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_pitchforks();
}
