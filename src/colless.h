#ifndef colless_h
#define colless_h


class colless_stat {

public:
  colless_stat(const std::vector< std::array< size_t, 2 >> & e) : edge(e) {
    tiplist = std::vector<int>(edge.size() + 2, -1);
  }

  size_t calc_colless() {
    size_t root_label = edge[0][0];

    std::sort(edge.begin(), edge.end(), [&](const auto& a, const auto& b) {
      return a[0] < b[0];
    });

    std::vector< size_t > s(edge.size());
    for (size_t i = 0; i < edge.size(); i++) {
      s[i] = get_diff_l_r(edge[i][1], root_label);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
  }

private:

  size_t get_diff_l_r(size_t label, size_t root_label) {
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

    std::vector< size_t > matches;
    auto match1 = std::lower_bound(edge.begin(), edge.end(), label, [&](const auto& a, size_t val){
      return a[0] < val;
    });

    while (match1 != edge.end()) {
      if ((*match1)[0] == label) {
        matches.push_back((*match1)[1]);
        match1++;
      } else {
        break;
      }
    }

    assert(matches.size() == 2);

    int l = get_num_tips(matches[0], root_label);
    int r = get_num_tips(matches[1], root_label);
    tiplist[matches[0]] = l;
    tiplist[matches[1]] = r;

    return std::abs(l - r);
  }

  size_t get_num_tips(size_t label, size_t root_label) {
    if (label >= tiplist.size()) {
      throw std::out_of_range("label > tiplist.size()");
    }

    if (label < root_label) {
      return 1;
    }

    if (tiplist[label] > 0) {
      return(tiplist[label]);
    }

    std::vector< size_t > matches;
    auto match1 = std::lower_bound(edge.begin(), edge.end(), label, [&](const auto& a, size_t val){
      return a[0] < val;
    });

    while (match1 != edge.end()) {
      if ((*match1)[0] == label) {
        matches.push_back((*match1)[1]);
        match1++;
      } else {
        break;
      }
    }

    if (matches.empty()) {
      // this can't really happen.
      tiplist[label] = 1;
      return 1;
    }

    size_t s = 0;
    for (size_t j = 0; j < matches.size(); ++j) {
      s += get_num_tips(matches[j], root_label);
    }
    tiplist[label] = s;
    return s;
  }

  std::vector< std::array< size_t, 2 >>  edge;
  std::vector<int> tiplist;
};



#endif
