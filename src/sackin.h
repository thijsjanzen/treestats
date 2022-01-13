#ifndef sackin_h
#define sackin_h

#include <vector>
#include <algorithm>

using ltable = std::vector< std::array<double, 4>>;


size_t calc_sackin(const ltable& ltable_) {
  std::vector< int > s_values(ltable_.size(), 0);
  s_values[0] = 1;
  s_values[1] = 1;

  std::vector< std::array<int, 2> > parent_map;

  for (int i = 0; i < ltable_.size(); ++i) {
    parent_map.push_back( {static_cast<int>(ltable_[i][2]), i});
  }

  std::sort(parent_map.begin(), parent_map.end(),
            [&](const auto& a, const auto& b) {
              return a[0] < b[0];
            });

  for (size_t i = 2; i < ltable_.size(); ++i) {
    int parent_id = ltable_[i][1];

    auto it = std::lower_bound(parent_map.begin(), parent_map.end(), parent_id,
                               [&](const auto& a, int ref) {
                                 return a[0] < ref;
                               });
    auto parent_index = (*it)[1];

    s_values[parent_index]++;
    s_values[i] = s_values[parent_index];
  }

  // verified with R for correct values
  // compared with apTreeShape results, based on ltable
  // 24-09-2021
  return(std::accumulate(s_values.begin(), s_values.end(), 0));
}

double correct_pda(const ltable& ltable_,
                  double Is) {
  double n = ltable_.size();
  double denom = powf(n, 1.5f);
  return Is / denom;;
}

double correct_yule(const ltable& ltable_,
                   double Is) {
  double n = ltable_.size();
  double sum_count = 0.f;
  for (double j = 2.f; j <= n; ++j) {
    sum_count += 1.f / j;
  }
  return (Is - 2 * n * sum_count) / n;
}

#endif
