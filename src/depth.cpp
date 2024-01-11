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
// DEPTH BASED statistics

#include <vector>
#include <array>
#include <Rcpp.h>

#include "util.h"        // NOLINT [build/include_subdir]
#include "depth.h"       // NOLINT [build/include_subdir]
#include "ltable.h"      // NOLINT [build/include_subdir]
#include "max_depth.h"   // NOLINT [build/include_subdir]
#include "b1.h"          // NOLINT [build/include_subdir]
#include "sym_nodes.h"   // NOLINT [build/include_subdir]

// [[Rcpp::export]]
int calc_max_del_width_cpp(const std::vector<int>& parent_list) {
  width::width_tree tree(parent_list);
  return tree.calc_max_del_width();
}

// [[Rcpp::export]]
double calc_max_del_width_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.max_del_width();
}

// [[Rcpp::export]]
int calc_max_width_cpp(const std::vector<int>& parent_list) {
  width::width_tree tree(parent_list);
  return tree.calc_max_width();
}

// [[Rcpp::export]]
double calc_max_width_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.calc_max_width();
}

// [[Rcpp::export]]
int calc_max_depth_cpp(const std::vector<int>& parent_list) {
  max_depth::max_depth_tree local_tree(parent_list);
  return local_tree.max_depth();
}

// [[Rcpp::export]]
double calc_max_depth_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.calc_max_depth();
}


// [[Rcpp::export]]
double calc_var_leaf_depth_cpp(const std::vector<int>& parent_list) {
  width::width_tree local_tree(parent_list);
  return local_tree.var_leaf_depth();
}

// [[Rcpp::export]]
double calc_var_leaf_depth_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat local_tree(l_in_cpp);
  return local_tree.calc_var_leaf_depth();
}

// [[Rcpp::export]]
int calc_sym_nodes_cpp(const std::vector<int>& parent_list) {
  sym_nodes::sym_node_tree focal_tree(parent_list);
  return focal_tree.calc_sym_nodes();
}

// [[Rcpp::export]]
double calc_b1_cpp(const std::vector<int>& parent_list) {
  b1_tree::b1_tree focal_tree(parent_list);
  return focal_tree.calc_b1();
}

// [[Rcpp::export]]
double calc_b1_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.calc_b1();
}

// [[Rcpp::export]]
double calc_b2_cpp(const std::vector<int>& parent_list) {
  width::width_tree tree(parent_list);
  return tree.calc_b2();
}

// [[Rcpp::export]]
double calc_b2_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.calc_b2();
}
