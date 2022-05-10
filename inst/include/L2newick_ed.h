#ifndef l_to_newick_ed_h
#define l_to_newick_ed_h

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <L2newick.h>

std::string ltable_to_newick_ed(const std::vector< std::array< double, 4>>& ltable, const double t,
                       bool drop_extinct) {
  auto L = ltable;
  // first sort ltable
  //  L = L[order(abs(L[, 3])), 1:4]
  std::sort(L.begin(), L.end(), [&](const auto& a, const auto& b) {
    return std::abs(static_cast<int>(a[2])) < std::abs(static_cast<int>(b[2]));
  });

  double age = t; // age = t
  std::vector< std::array< double, 4>> new_L;

  for (auto& i : L) {
    // test i[3] == -1 but avoid precision issue
    bool is_extant = ((i[3] + 1) < 0.000001);
    if (is_extant) {
      i[3] = age;
      if (drop_extinct == true) {
        new_L.push_back(i);
      }
    }
  }
  
  // keep a copy of the original ltable for later lookup purpose
  auto L_original = L;
  L_original[0][0] = -1.0;

  if (drop_extinct == true) {
    L = new_L;
  } else {
    // L[0][0] cannot be -1 when extinct lineages are dropped
    L[0][0] = -1.0;
  }

  std::vector< std::string > linlist_4(L.size());
  size_t index = 0;
  for (const auto& i : L) {
    std::string add = "t" + std::to_string(abs(static_cast<int>(i[2])));
    linlist_4[index] = add;
    index++;
  }

  // verified correct for 4/5 tip phylogeny up until here.
  while(true) {
    auto j = which_max_index(L);
    int parent    = static_cast<int>(L[j][1]);
    int parentj   = index_of_parent(L, parent);
    if (parentj != -1) {
      double bl = std::abs(L[parentj][3] - L[j][0]);
      std::string spec1 = linlist_4[parentj] + ":" + d_to_s(bl);
      double bl2 = std::abs(L[j][3] - L[j][0]);
      std::string spec2 = linlist_4[j] + ":" + d_to_s(bl2);
      linlist_4[parentj] = "(" + spec1 + "," + spec2 + ")";
      L[parentj][3] = L[j][0];
      remove_from_dataset(L, linlist_4, j);
    } else {
      parentj = index_of_parent(L_original, parent);
      for (auto i = 0; i <= 2; ++i) {
        L[j][i] = L_original[parentj][i];
      }
    }

    if (linlist_4.size() == 1) {
      break;
    }
  }

  return linlist_4[0] + ":" + d_to_s(L[0][3]) + ";";
}




#endif
