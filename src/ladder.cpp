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

#include <vector>
#include <array>
#include <RcppArmadillo.h>

#include "util.h"        // NOLINT [build/include_subdir]
#include "avgladder.h"   // NOLINT [build/include_subdir]

// [[Rcpp::export]]
double avgLadder_cpp(const std::vector<int>& tree_edge) {
  try {
    return calc_ladder(tree_edge, false);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// [[Rcpp::export]]
double max_ladder_cpp(const std::vector<int>& tree_edge) {
  try {
    return calc_ladder(tree_edge, true);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

