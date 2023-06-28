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
#pragma once

#include <vector>

inline double binom_coeff(const int& n, const int& k) {
  if (n == k) return(0);

  if (n < 100) {
    std::vector<int> aSolutions(k);
    aSolutions[0] = n - k + 1;

    for (int i = 1; i < k; ++i) {
      aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
    }

    return aSolutions[k - 1];
  } else {
    double s = 0.0;
    for (size_t i = 1; i < k; ++i) {
      s += std::log(static_cast<double>(n - i + 1));
      s -= log(i);
    }
    return(exp(s));
  }
}

inline double log_binom_coeff(const int n, const int k) {
  double s = 0.0;
  for (size_t i = 1; i < k; ++i) {
    s += std::log(static_cast<double>(n - i + 1));
    s -= log(i);
  }
  return(s);
}




inline double binom_coeff_2(int n) {
  double n_d = static_cast<double>(n);
  return (n_d - 1) * n_d * 0.5;
}
