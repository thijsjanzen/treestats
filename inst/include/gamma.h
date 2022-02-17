#ifndef gamma_h
#define gamma_h

#include <vector>
#include <numeric>

const double prefactor = 2 * sqrt(3);

double calc_gamma(const std::vector<double>& brts_in) {
  std::vector<double> brts_ = brts_in;
  double n = brts_in.size() + 1;

  auto h = brts_[0]; //*std::max_element(brts_.begin(), brts_.end());

  for (auto& i : brts_) {
    i =  h - i;
  }

  std::sort(brts_.begin(), brts_.end());

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


  double mult_total = 1.0 / total;
  return prefactor * sqrt(n - 2) * (double_sum * 1.0 / (n - 2) - total * 0.5) * mult_total;
}

#endif
