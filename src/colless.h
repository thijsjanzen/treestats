#ifndef colless_h
#define colless_h

#include <vector>
#include <array>
#include <numeric> // std::accumulate

class colless_stat {

public:
  colless_stat(const std::vector< int>& p,
                size_t n_tips) : parents(p), num_tips(n_tips) {
  }

  size_t calc_colless() {
    tiplist = std::vector< int >(parents.size(), 0);
    for (size_t i = 1; i <= num_tips; ++i) {
      tiplist[ parents[i] ]++;
    }

    size_t s = 0;
    for (size_t i = tiplist.size() - 1; i > num_tips + 1; i--) {
      if (tiplist[ parents[i] ] > 0) {
        int l = tiplist[ parents[i] ];
        int r = tiplist[i];
        l - r < 0 ? s -= l - r : s+= l - r;
      }
      tiplist[ parents[i] ] += tiplist[i];
    }
    return s;
  }

  double correct_pda(double Ic) {
    double denom = powf(num_tips, 1.5f);
    return 1.0 * Ic / denom;
  }

  double correct_yule(double Ic) {
    static const double g = 0.577215664901532;
    auto output = (Ic - num_tips * log(num_tips) - num_tips * (g - 1 - log(2))) / num_tips;
    return output;
  }


private:
  const std::vector< int >  parents;
  std::vector< int > tiplist;
  const size_t num_tips;
};



#endif
