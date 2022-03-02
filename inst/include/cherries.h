#ifndef CHERRIES_H
#define CHERRIES_H

#include <vector>
#include <array>

using ltable = std::vector< std::array<double, 4>>;

size_t find_daughters(const ltable& ltab_in,
                      double id,
                      double bt) {
  size_t num_daughters = 0;

  auto start = std::lower_bound(ltab_in.begin(), ltab_in.end(), bt,
                                [](const auto& a, double val) {
                                  return a[0] > val; });

  for (auto i = start; i != ltab_in.end(); ++i) {
    if ((*i)[1] == id && (*i)[0] <= bt) {
      num_daughters++;
      if (num_daughters > 1) break;
    }
  }
  return num_daughters;
}

size_t calc_cherries_ltable(const ltable& ltab_in) {
  size_t num_cherries = 0;
  for (const auto& i : ltab_in) {
 //   if (i[3] != -1) continue;  // non-extant species can not be a cherry

    auto parent = i[1];

    if (parent == 0) continue;
  //  if (ltab_in[std::abs(parent)][3] != -1) continue;
    auto bt = i[0];

    size_t num_daughter_branches = find_daughters(ltab_in, parent, bt);
    size_t own_daughters         = find_daughters(ltab_in, i[2], bt);

    if (num_daughter_branches == 1 && own_daughters == 0) {
      num_cherries++;
    }
  }
  return num_cherries;
}



#endif



