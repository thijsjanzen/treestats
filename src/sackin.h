#ifndef sackin_h
#define sackin_h

#include <vector>

using ltable = std::vector< std::array<double, 4>>;

template <class IT>
size_t find_parent(IT first,
                   IT last,
                   int focal_id,
                   int start_index) {

  auto focal = first + start_index;
  for (; focal >= first; --focal) {
    if ( (*focal)[1] == focal_id )
      return std::distance(first, focal);
  }

  if (focal != first) {
    // force_output("could not find parent, retrying\n");
    return find_parent(first, last, focal_id, std::distance(first, last));
  } else {
    //  force_output("could not find parent at all\n");
    throw std::out_of_range("could not find parent\n");
  }
  return -1;
}


size_t calc_sackin(const ltable& ltable_) {
  std::vector< int > s_values(ltable_.size(), 0);
  s_values[0] = 1;
  s_values[1] = 1;

  std::vector< std::array<int, 2>> ref_tab(ltable_.size());
  for (size_t i = 0; i < ltable_.size(); ++i) {
    ref_tab[i] = {static_cast<int>(ltable_[i][1]),
                  static_cast<int>(ltable_[i][2])};
  }

  for (size_t i = 2; i < ref_tab.size(); ++i) {
    auto parent_id = ref_tab[i][0];
    auto parent_index = find_parent(ref_tab.begin(), ref_tab.end(), parent_id, i);
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
