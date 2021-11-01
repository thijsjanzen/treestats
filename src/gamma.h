#ifndef gamma_h
#define gamma_h

#include <vector>

struct gamma_stat {

  gamma_stat(const std::vector<float>& brts_in) : brts_(brts_in), n(brts_in.size() + 1) {
    auto h = *std::max_element(brts_.begin(), brts_.end());
    for (auto& i : brts_) {
      i =  h - i;
    }
    brts_.push_back(0.f);
    brts_.push_back(h);
    std::sort(brts_.begin(), brts_.end());

    std::adjacent_difference(brts_.cbegin(), brts_.cend(), std::back_inserter(g));
  }

  double calc_gamma_stat() {

    float double_sum = 0.f;
    float total = 0.f;
    for (size_t i = 1; i <= n; ++i) {

      if (i >= 2) {
        total += i * g[i];
      }

      if (i <= (n-1)) {
        for (size_t k = 1; k <= i; ++k) {
          double_sum += k * g[k];
        }
      }
    }

    float a = double_sum * 1.0f / (n - 2);

    float b = total / 2;
    float c = total * sqrtf(1.f / (12 * n - 24));

    return (a - b) / c;
  }



private:
  std::vector<float> brts_;
  const float n;
  std::vector<float> g;

};

#endif
