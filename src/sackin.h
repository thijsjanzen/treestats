#ifndef sackin_h
#define sackin_h

#include <vector>

using ltable = std::vector< std::array<double, 4>>;

size_t find_parent(const ltable& ltable_,
                   int focal_id,
                   int start_index) {

  for(int i = start_index; i >= 0; i--) {
    if (static_cast<int>(ltable_[i][2]) == focal_id) {
      return i;
    }
  }

  if (start_index != ltable_.size()) {
   // force_output("could not find parent, retrying\n");
    return find_parent(ltable_, focal_id, ltable_.size());
  } else {
  //  force_output("could not find parent at all\n");
    return -1; // trigger access violation --> update to throw
  }
}


size_t calc_sackin(const ltable& ltable_) {
  std::vector< int > s_values(ltable_.size(), 0);
  s_values[0] = 1;
  s_values[1] = 1;

  // ltable:
  // 0 = branching time // not used here
  // 1 = parent
  // 2 = id
  // 3 = extinct time // not used here
  for (size_t i = 2; i < ltable_.size(); ++i) {
    int parent_id = ltable_[i][1];
    int parent_index = find_parent(ltable_, parent_id, i);
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
