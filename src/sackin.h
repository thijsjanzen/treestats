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

    size_t s = 0;
    for (int i = edge.size() - 1; i >= 0; i--) {
      s += get_num_tips(edge[i][1], root_label);
    }

    return s;
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


class sackin_stat2 {

public:
  sackin_stat2(const std::vector< int>& p,
               size_t n_tips) : parents(p), num_tips(n_tips) {
  }



  size_t calc_sackin() {
    tiplist = std::vector< int >(parents.size(), 0);
    for (size_t i = 1; i <= num_tips; ++i) {
      tiplist[ parents[i]]++;
    }

    for (size_t i = tiplist.size() - 1; i > num_tips + 1; i--) {

      tiplist[ parents[i] ] += tiplist[i];
    }

    return std::accumulate(tiplist.begin(), tiplist.end(), 0);
  }

  double calc_blum() {
    tiplist = std::vector< int >(parents.size(), 0);
    for (size_t i = 1; i <= num_tips; ++i) {
      tiplist[ parents[i]]++;
    }

    for (size_t i = tiplist.size() - 1; i > num_tips + 1; i--) {

      tiplist[ parents[i] ] += tiplist[i];
    }

    double s = 0.0;
    for (size_t i = num_tips + 1; i < tiplist.size(); ++i) {
      s += log(1.0 * tiplist[i]);
    }
    return s;
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
  const std::vector< int >  parents;
  std::vector< int > tiplist;
  const size_t num_tips;
};

using ltable = std::vector< std::array<double, 4>>;

class sackin_stat_ltab {
public:
  sackin_stat_ltab(const ltable& ltab_in) : ltable_(ltab_in) {
  }

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

  size_t calc_sackin() {
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

  double correct_pda(double Is) {
    size_t n = ltable_.size();
    double denom = powf(n, 1.5f);
    return 1.0 * Is / denom;
  }

  double correct_yule(double Is) {
    double sum_count = 0.0;
    size_t n = ltable_.size();
    for (size_t j = 2; j <= n; ++j) {
      sum_count += 1.0 / j;
    }
    return 1.0 * (Is - 2.0 * n * sum_count) / n;
  }

private:
  const ltable ltable_;
};




#endif
