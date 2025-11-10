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
//
#include <vector>
#include <array>
#include <string>
#include <RcppArmadillo.h>

#include "util.h"        // NOLINT [build/include_subdir]
#include "beta.h"        // NOLINT [build/include_subdir]


// [[Rcpp::export]]
double calc_beta_cpp(const Rcpp::List& phy,
                     double upper_lim,
                     std::string algorithm,
                     double abs_tol,
                     double rel_tol) {
  try {
    Rcpp::NumericMatrix edge = phy["edge"];
    if (edge.nrow() == 2) {
      Rcpp::warning("Trees with only two tips have undefined beta");
      return NA_REAL;
    }

    std::vector< std::array< int, 2 >> local_edge(edge.nrow());
    for (int i = 0; i < edge.nrow(); ++i) {
      local_edge[i] = {static_cast<int>(edge(i, 0)),
                       static_cast<int>(edge(i, 1))};
    }

    return calc_beta(local_edge, -2.0, upper_lim, algorithm, abs_tol, rel_tol);
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
double calc_beta_ltable_cpp(const Rcpp::NumericMatrix& ltable,
                            double upper_lim,
                            std::string algorithm,
                            double abs_tol,
                            double rel_tol) {
  try {
    std::vector< std::array< double, 4 >> ltab(ltable.nrow());
    for (int i = 0; i < ltable.nrow(); ++i) {
      std::array< double, 4> row_entry = {ltable(i, 0), ltable(i, 1),
                                          ltable(i, 2), ltable(i, 3)};
      ltab[i] = row_entry;
    }

    return calc_beta(ltab, -2.0, upper_lim, algorithm, abs_tol, rel_tol);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
