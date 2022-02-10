#ifndef colless_h
#define colless_h

#include <vector>
#include <array>
#include <numeric> // std::accumulate

#include <iostream>


using ltable = std::vector< std::array<double, 4>>;

class colless_stat_ltable {
public:
  colless_stat_ltable(const ltable& l_in) : ltable_(l_in) {
    extant_tips = std::vector<int>(l_in.size(), 1);
    num_tips = get_num_tips();
  }

  size_t calc_colless() {
    size_t colless_stat = 0;
    while(true) {
      auto j = get_min_index();
      auto parent = ltable_[j][1];
      if (parent == 0) {// we hit the root!
        j++;
        parent = ltable_[j][1];
      }
      auto j_parent = index_of_parent(parent);

      int L = extant_tips[j];
      int R = extant_tips[j_parent];
      colless_stat += std::abs(L - R);
      extant_tips[j_parent] = L + R;
      remove_from_dataset(j);

      if (ltable_.size() == 1) break;
    }
    return colless_stat;
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

  size_t get_min_index() {
    auto min_val = std::min_element(ltable_.begin(), ltable_.end(),
                                    [&](const auto& a, const auto& b) {
                                      return a[0] < b[0];
                                    });
    return std::distance(ltable_.begin(), min_val);
  }

  int index_of_parent(int parent) {
    int index = 0;
    bool found = false;
    for (; index < ltable_.size(); ++index) {
      if (ltable_[index][2] == parent) {
        found = true;
        break;
      }
    }
    if (!found) index = -1;
    return index;
  }

  void remove_from_dataset(size_t index) {
    std::swap(extant_tips[index], extant_tips.back());
    extant_tips.pop_back();
    std::swap(ltable_[index], ltable_.back());
    ltable_.pop_back();
  }

  size_t get_num_tips() {
   /* size_t num_extant_tips = 0;
    for (const auto& i : ltable_) {
      if (i[3] < 0) num_extant_tips++;
    }
    return num_extant_tips;*/
   return ltable_.size();
  }
  ltable ltable_;
  std::vector< int > extant_tips;
  size_t num_tips;
};



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
