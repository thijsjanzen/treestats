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

// we don't lint this file, it generates many weird lints on tab size / number
// of spaces
//
// NOLINTBEGIN

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector Ax_tree(const Rcpp::IntegerMatrix &edge,
                            const Rcpp::NumericVector &lengths,
                            const Rcpp::NumericVector &x,
                            int nNodes) {
   Rcpp::NumericVector result(nNodes);
   int nEdges = edge.nrow();
   for (int i = 0; i < nEdges; i++) {
      int u = edge(i, 0) - 1;            // convert to zero-based
      int v = edge(i, 1) - 1;
      double w = lengths[i];
      // adjacency contribution
      result[u] += w * x[v];
      result[v] += w * x[u];
   }
   return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector Lx_tree_weighted(const Rcpp::IntegerMatrix &edge,
                                     const Rcpp::NumericVector &w,
                                     const Rcpp::NumericVector &x,
                                     int nNodes) {
   Rcpp::NumericVector result(nNodes);
   int nEdges = edge.nrow();
   for (int i = 0; i < nEdges; i++) {
      int u = edge(i, 0) - 1;  // zero-based
      int v = edge(i, 1) - 1;
      double weight = w[i];
      double diff_uv = x[u] - x[v];
      result[u] += weight * diff_uv;
      result[v] -= weight * diff_uv;  // = weight * (x[v] - x[u])
   }
   return result;
}

// NOLINTEND
