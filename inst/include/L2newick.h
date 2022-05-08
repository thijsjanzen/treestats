#ifndef l_to_newick_h
#define l_to_newick_h

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>


size_t which_max_index(const std::vector< std::array< double, 4>>& ltable) {
  auto max_val = std::max_element(ltable.begin(), ltable.end(),
                                  [&](const auto& a, const auto& b) {
                                    return a[0] < b[0];});
  return std::distance(ltable.begin(), max_val);
}

int index_of_parent(const std::vector< std::array< double, 4>>& ltable,
                    int parent) {
  int index = 0;
  bool found = false;
  for (; index < ltable.size(); ++index) {
    if (ltable[index][2] == parent) {
      found = true;
      break;
    }
  }
  if (!found) index = -1;
  return index;
}

std::string d_to_s(double d) {
  std::stringstream out;
  out << std::fixed << std::setprecision(15) << d;
  return out.str();
}

void remove_from_dataset(std::vector< std::array< double, 4>>& ltable,
                         std::vector< std::string>& linlist,
                         size_t index) {
  std::swap(ltable[index], ltable.back());
  ltable.pop_back();
  std::swap(linlist[index], linlist.back());
  linlist.pop_back();
}

std::string ltable_to_newick(const std::vector< std::array< double, 4>>& ltable,
                       bool drop_extinct) {
  auto L = ltable;
  // first sort ltable
  //  L = L[order(abs(L[, 3])), 1:4]
  std::sort(L.begin(), L.end(), [&](const auto& a, const auto& b) {
    return std::abs(static_cast<int>(a[2])) < std::abs(static_cast<int>(b[2]));
  });

  double age = L[0][0]; // age = L[1, 1]
  std::vector< std::array< double, 4>> new_L;

  for (auto& i : L) {
    i[0] = age - i[0]; // L[, 1] = age - L[, 1]

    bool is_extant = i[3] == -1;

    if (i[3] != -1) { // notmin1 = which(L[, 4] != -1)
      // L[notmin1, 4] = age - L[notmin1, 4]
      i[3] = age - i[3];
    } else {
      i[3] = age;
    }

    if (drop_extinct == true) {
      if (is_extant) {
        new_L.push_back(i);
      }
    }
  }

  if (drop_extinct == true) {
    L = new_L;
  }

  L[0][0] = -1.0;

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
      double bl = L[parentj][3] - L[j][0];
      std::string spec1 = linlist_4[parentj] + ":" + d_to_s(bl);
      double bl2 = L[j][3] - L[j][0];
      std::string spec2 = linlist_4[j] + ":" + d_to_s(bl2);
      linlist_4[parentj] = "(" + spec1 + "," + spec2 + ")";
      L[parentj][3] = L[j][0];
      std::swap(linlist[j], linlist.back());
      linlist.pop_back();
    } else {
      for (auto i = 0; i <= 2; ++i) {
        L[j][i] = L[parentj][i];
      }
    }

    if (linlist_4.size() == 1) {
      break;
    }
  }

  return linlist_4[0] + ":" + d_to_s(L[0][3]) + ";";
}




#endif
