#ifndef gamma_h
#define gamma_h

#include <vector>
double calc_gamma(const std::vector<double>& brts_in) {
  std::vector<double> brts_ = brts_in;

  auto h = *std::max_element(brts_.begin(), brts_.end());

  for (auto& i : brts_) {
    i =  h - i;
  }

  brts_.push_back(0.0);
  brts_.push_back(h);
  std::sort(brts_.begin(), brts_.end());

  double total = 0.0;
  double g = brts_[1] - brts_[0];
  double local_double_sum = g;
  double double_sum = g;

  for (size_t i = 2; i < brts_.size(); ++i) {
    g = brts_[i] - brts_[i - 1];

    total += i * g;

    local_double_sum += i * g;
    double_sum += local_double_sum;
  }

  double n = brts_in.size() + 1;
  if (n > brts_.size()) {
    throw std::out_of_range("n > brts_.size()");
  }
  g = brts_[n] - brts_[n - 1];
  total += n * g;


  double a = double_sum * 1.0 / (n - 2);

  double b = total / 2;
  double c = total * sqrt(1.f / (12 * n - 24));

  return (a - b) / c;
}



#endif
