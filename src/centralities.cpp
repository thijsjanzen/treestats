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
// Network statistics based

#include <vector>
#include <array>
#include <RcppArmadillo.h>

#include "centralities.h"  // NOLINT [build/include_subdir]
#include "util.h"          // NOLINT [build/include_subdir]
#include "dist_nodes.h"    // NOLINT [build/include_subdir]

// [[Rcpp::export]]
double calc_wiener_cpp(const std::vector<int>& edge_in,
                       const std::vector<double>& el,
                       bool normalize,
                       bool weight) {
  edge e = flat_to_table(edge_in);
  return wiener(e, el, normalize, weight, is_binary(edge_in));
 // return is_binary(edge_in) ?
   //       wiener(e, el, normalize, weight)  :
  //        wiener_poly(e, el, normalize, weight);
}

// [[Rcpp::export]]
double calc_max_betweenness_cpp(const std::vector<int>& edge_in,
                                const std::vector<double>& el) {
  edge e = flat_to_table(edge_in);

 // return is_binary(edge_in) ? max_betweenness(e, el) :
  //                            max_betweenness_poly(e, el);
  return max_betweenness(e, el, is_binary(edge_in));
}

// [[Rcpp::export]]
double calc_max_betweenness_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  return max_betweenness_ltable(l_in_cpp);
}

// [[Rcpp::export]]
double calc_max_closeness_cpp(const std::vector<int>& edge_in,
                              const std::vector<double>& el,
                              bool weight) {

  edge e = flat_to_table(edge_in);

 // return is_binary(edge_in) ?
//             1.0 / min_farness(e, el, weight) :
//             1.0 / min_farness_poly(e, el, weight);
  return 1.0 / min_farness(e, el, weight, is_binary(edge_in));
}

// [[Rcpp::export]]
double calc_diameter_cpp(const std::vector<int>& edge_in,
                         const std::vector<double>& el,
                         bool weight) {
  edge e = flat_to_table(edge_in);

 // return is_binary(edge_in) ?
 //               diameter(e, el, weight) :
  //              diameter_poly(e, el, weight);
  return diameter(e, el, weight, is_binary(edge_in));
}

// [[Rcpp::export]]
double calc_diameter_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                                bool weight) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  return diameter_ltable(l_in_cpp, weight);
}
