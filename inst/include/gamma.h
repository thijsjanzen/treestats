#ifndef gamma_h
#define gamma_h

#include <vector>
#include <numeric>

double calc_gamma(std::vector<double> brts_) {
  double n = brts_.size() + 1;

  auto h = brts_[0]; // assuming brts are in T...0

  for (auto& i : brts_) {
    i =  h - i;
  }

  double total = 0.0;
  double double_sum = 0.0;
  double temp = n * (h - brts_.back());
  std::adjacent_difference(brts_.begin(), brts_.end(), brts_.begin());

  size_t j = 1;
  for(auto i : brts_) {
    total += j * i;
    double_sum += total;
    j++;
  }

  total += temp;

  const double prefactor = 2 * sqrt(3);
  double mult_total = 1.0 / total;
  return prefactor * sqrt(n - 2) * (double_sum * 1.0 / (n - 2) - total * 0.5) * mult_total;
}

#endif
