#ifndef binom_h
#define binom_h

#include <vector>

inline int binom_coeff(const int& n, const int& k) {
  std::vector<int> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (int i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

inline int binom_coeff_2(int n) {
   return (n - 1) * n * 0.5;
}

#endif
