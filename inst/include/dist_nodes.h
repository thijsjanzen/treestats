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
#pragma once

#include <vector>
#include <array>
#include <Rcpp.h>

// The function below is lifted from the ape package.
inline std::vector< std::vector< double >> dist_nodes(
    const std::vector< std::array< size_t, 2 >>& edge,
    const std::vector<double>& el) {

  int n = 1 + edge.size() / 2;
  int m = n - 1;
  auto nm = n + m;
  static double max_s = 46340;   // floor(sqrt(2^31 - 1))
  if (nm > max_s) {
     throw std::runtime_error("tree too big");
  }
  // code below is from the Ape package
  std::vector< size_t > e1(edge.size());
  std::vector< size_t > e2(edge.size());

  for (size_t i = 0; i < edge.size(); ++i) {
    e1[i] = edge[i][0] - 1;
    e2[i] = edge[i][1] - 1;
  }

  int k = 0, a, d, NM = n + m, ROOT;
  double x;
  size_t N = e1.size();
  std::vector< std::vector<double>> D(NM, std::vector<double>(NM, 0.0));

  ROOT = e1[0]; d = e2[0]; /* the 2 nodes of the 1st edge */
  D[ROOT][d] = D[d][ROOT] = el[0]; /* the 1st edge gives the 1st distance */

  /* go down along the edge matrix
   starting at the 2nd edge: */
  for (size_t i = 1; i < N; i++) {
    a = e1[i]; d = e2[i]; x = el[i]; /* get the i-th nodes and branch length */
  D[a][d] = D[d][a] = x;
  /* then go up along the edge matrix from the i-th edge
   to visit the nodes already visited and update the distances: */
  for (int j = i - 1; j >= 0; j--) {
    k = e2[j];
    if (k == a) continue;
    D[k][d] = D[d][k] = D[a][k] + x;
  }
  if (k != ROOT)
    D[ROOT][d] = D[d][ROOT] = D[ROOT][a] + x;
  }
  return D;
}
