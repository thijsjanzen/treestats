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
//
// SACKIN BASED Summary statistics

#include <vector>
#include <array>
#include <Rcpp.h>

#include "util.h"     // NOLINT [build/include_subdir]
#include "sackin.h"   // NOLINT [build/include_subdir]
#include "ltable.h"   // NOLINT [build/include_subdir]
#include "cherries.h" // NOLINT [build/include_subdir]

// [[Rcpp::export]]
double calc_sackin_cpp(const std::vector<int>& tree_edge,
                       const Rcpp::String& normalization) {
  sackin::sackin_tree sackin_tree(tree_edge);
  double output = static_cast<double>(sackin_tree.calc_sackin());

  if (normalization == "yule") {
    size_t n  = tree_edge.size() / 4 + 1;
    output = correction::correct_yule(n, output);
  }
  if (normalization == "pda") {
    size_t n  = tree_edge.size() / 4 + 1;
    output = correction::correct_pda(n, output);
  }

  return output;
}


// [[Rcpp::export]]
double calc_sackin_ltable_cpp(const Rcpp::NumericMatrix& ltab,
                       const Rcpp::String& normalization) {
  auto local_ltab = convert_to_ltable(ltab);
  return calc_sackin(local_ltab, normalization);
}


// [[Rcpp::export]]
double calc_tot_coph_cpp(const std::vector<int>& tree_edge) {
  sackin::sackin_tree sackin_tree(tree_edge);
  return sackin_tree.calc_tot_coph();
}

// [[Rcpp::export]]
double calc_tot_coph_ltable_cpp(const Rcpp::NumericMatrix& ltab) {
  auto local_ltab = convert_to_ltable(ltab);
  ltab::stat s(local_ltab);
  return s.calc_tot_coph();
}

// [[Rcpp::export]]
double calc_blum_cpp(const std::vector<int>& tree_edge,
                     bool normalize) {
  sackin::sackin_tree sackin_tree(tree_edge);
  double output = sackin_tree.calc_blum();
  if (normalize)  {
    size_t n  = tree_edge.size() / 4 + 1;
    output = correction::correct_blum(n, output);
  }
  return output;
}
// [[Rcpp::export]]
double calc_blum_ltable_cpp(const Rcpp::NumericMatrix& ltab_in,
                            bool normalize) {
  auto local_ltab = convert_to_ltable(ltab_in);
  ltab::stat s(local_ltab);
  return s.calc_blum(normalize);
}


// [[Rcpp::export]]
size_t cherries_cpp(const std::vector<int>& tree_edge) {
  sackin::sackin_tree sackin_tree(tree_edge);
  return sackin_tree.count_cherries();
}

// [[Rcpp::export]]
size_t cherries_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  return calc_cherries_ltable(local_ltab);
}

// [[Rcpp::export]]
size_t pitchforks_cpp(const std::vector<int>& tree_edge) {
  // ltable version uses colless
  sackin::sackin_tree sackin_tree(tree_edge);
  return sackin_tree.count_pitchforks();
}

