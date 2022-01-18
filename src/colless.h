#ifndef colless_h
#define colless_h

#include <vector>
#include <array>
#include <numeric> // std::accumulate


class colless_stat {

public:
  colless_stat(const std::vector< std::array< size_t, 2 >> & e) : edge(e) {
    tiplist = std::vector<size_t>(edge.size() + 2, 0);
  }

  size_t calc_colless() {
    size_t root_label = edge[0][0];

    std::sort(edge.begin(), edge.end(), [&](const auto& a, const auto& b) {
      return a[0] < b[0];
    });

    std::vector< size_t > s;
    for (size_t i = 0; i < (edge.size() - 1); i++) {

      get_num_tips(edge[i][1], root_label);
      if (i + 1 < edge.size()) {
        if (edge[i][0] == edge[i + 1][0]) {
          get_num_tips(edge[i + 1][1], root_label);

          int l = static_cast<int>(tiplist[ edge[i][1]    ]);
          int r = static_cast<int>(tiplist[ edge[i + 1][1]]);
          int dd = l - r;

          s.push_back(std::abs(dd));
          i++;
        }
      }
    }
    return  std::accumulate(s.begin(), s.end(), 0.0);
  }

  double correct_pda(size_t n,
                     double Ic) {
    double denom = powf(n, 1.5f);
    return 1.0 * Ic / denom;
  }

  double correct_yule(size_t n,
                      double Ic) {
    static const double g = 0.577215664901532;
    auto output = (Ic - n * log(n) - n * (g - 1 - log(2))) / n;
    return output;
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
  std::vector<size_t> tiplist;
};



#endif
