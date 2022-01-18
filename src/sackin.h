#ifndef sackin_h
#define sackin_h

#include <vector>
#include <algorithm>
#include <array>

#include "Rcpp.h"

class sackin_stat {

public:
  sackin_stat(const std::vector< std::array< size_t, 2 >> & e) : edge(e) {
    tiplist = std::vector<int>(edge.size() + 2, -1);
  }

  size_t calc_sackin() {
    size_t root_label = edge[0][0];

    std::sort(edge.begin(), edge.end(), [&](const auto& a, const auto& b) {
      return a[0] < b[0];
    });

    std::vector< size_t > s(edge.size());
    for (size_t i = 0; i < edge.size(); i++) {
      s[i] = get_num_tips(edge[i][1], root_label);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
  }


  double correct_pda(size_t n,
                     double Is) {
    double denom = powf(n, 1.5f);
    return 1.0 * Is / denom;
  }

  double correct_yule(size_t n,
                      double Is) {
    double sum_count = 0.0;
    for (size_t j = 2; j <= n; ++j) {
      sum_count += 1.0 / j;
    }
    return 1.0 * (Is - 2.0 * n * sum_count) / n;
  }

private:

  size_t get_num_tips(size_t label, size_t root_label) {
    if (label >= tiplist.size()) {
      throw std::out_of_range("label > tiplist.size()");
    }

    if (label < root_label) {
      tiplist[label] = 1;
      return 1;
    }

    if (tiplist[label] > 0) {
      return(tiplist[label]);
    }

    std::vector< size_t > matches(2);
    auto match1 = std::lower_bound(edge.begin(), edge.end(), label, [&](const auto& a, size_t val){
      return a[0] < val;
    });

    if (match1 != edge.end()) {
      if ((*match1)[0] == label) {
        matches[0] = (*match1)[1];
        match1++;
        if ((*match1)[0] == label) {
          matches[1] = (*match1)[1];
        } else {
          matches.pop_back();
        }
      }
    } else {
      // this can't really happen.
      tiplist[label] = 1;
      return 1;
    }

    size_t s = 0;
    for (auto i : matches) {
      s += get_num_tips(i, root_label);
    }
    tiplist[label] = s;
    return s;
  }

  std::vector< std::array< size_t, 2 >>  edge;
  std::vector<int> tiplist;


};





#endif
