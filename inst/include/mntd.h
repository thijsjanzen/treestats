#ifndef mntd_h
#define mntd_h

#include <vector>
#include <array>
#include <numeric>

using ltable = std::vector< std::array<double, 4>>;

struct lower_tri {
  lower_tri(size_t n) : n_(n) {
    data_ = std::vector<double>((n_ * (n_ - 1)) * 0.5, 0.0);
  }

  int get_linear_index(int i, int j) {
    int index = i > j ? i * (i - 1) * 0.5 + j :
                        j * (j - 1) * 0.5 + i;

    if (index < 0) index = 0;
    return index;
  }

  void set_val(int i, int j, double val) {
    if (i == j) return; // do nothing, these values don't exist.

    auto local_index = get_linear_index(i, j);

    if (local_index < 0 || local_index > data_.size()) {
      throw "local_index outside data_";
    }

    data_[local_index] = val;
  }

  double get_val(int i, int j) {
    if (i == j) return 0.0;

    auto local_index = get_linear_index(i, j);

    if (local_index < 0 || local_index > data_.size()) {
      throw "local_index outside data_";
    }
     return data_[local_index];
  }

  double get_sum_tips() {
    return std::accumulate(data_.begin(), data_.begin() + n_ - 1, 0.0);
  }

  std::vector<double> data_;
  size_t n_;
};

lower_tri dist_nodes_tri(const std::vector< std::array< size_t, 2 >>& edge,
                         const std::vector<double>& el) {
  int n = 1 + edge.size() / 2;
  int m = n - 1;
  auto nm = n + m;
  static double max_s = 46340; // floor(sqrt(2^31 - 1))
  if (nm > max_s) {
    throw std::runtime_error("tree too big");
  }
  // code below is from the Ape package
  int i, j, k, a, d, NM = n + m, ROOT;
  double x;
  size_t N = edge.size();
  lower_tri D(NM);

  ROOT = edge[0][0] - 1;
  d    = edge[0][1] - 1; /* the 2 nodes of the 1st edge */

  D.set_val(d, ROOT, el[0]);

  /* go down along the edge matrix
   starting at the 2nd edge: */
  for (i = 1; i < N; i++) {
      a = edge[i][0] - 1;
      d = edge[i][1] - 1;
      x = el[i]; /* get the i-th nodes and branch length */
      D.set_val(a, d, x);
      /* then go up along the edge matrix from the i-th edge
       to visit the nodes already visited and update the distances: */
      for (j = i - 1; j >= 0; j--) {
        k = edge[j][1] - 1;
        if (k == a) continue;

        double val_to_set = D.get_val(a, k) + x;
        D.set_val(k, d, val_to_set);
      }
      if (k != ROOT) {
        double val_to_set = D.get_val(ROOT, a) + x;
        D.set_val(ROOT, d, val_to_set);
      }
  }
  return D;
}

double calc_mntd_ltable(const ltable& ltable_) {
  std::vector<double> dist(ltable_.size() + 1, -1);

  for (const auto& i : ltable_) {
    auto parent = std::abs(i[1]);
    auto daughter = std::abs(i[2]);
    auto dist_to_nearest_taxon = 2 * i[0];
    if (i[3] != -1) {
      dist_to_nearest_taxon = i[0] + (i[0] - i[3]);
    }

    if (daughter < 0 || daughter > dist.size()) {
      throw std::out_of_range("daughter outside dist");
    }

    dist[daughter] = dist_to_nearest_taxon;

    if (parent < 0 || parent > dist.size()) {
      throw std::out_of_range("parent outside dist");
    }

    if (dist[parent] > 0) {
      if (dist_to_nearest_taxon < dist[parent]) {
        dist[parent] = dist_to_nearest_taxon;
      }
    } else {
      dist[parent] = dist_to_nearest_taxon;
    }
  }

  dist[0]= 0.0;
  auto sum_dist = std::accumulate(dist.begin(), dist.end(), 0.0);
  return sum_dist * 1.0 / ltable_.size();
}

double calc_mntd_stat(const std::vector< std::array< size_t, 2 >>& edge,
                      const std::vector<double>& el) {
  auto dist_mat = dist_nodes_tri(edge, el);

  int n = (el.size() + 2) * 0.5;

  double mntd = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double min_val = 1e6;
    for (size_t j = 0; j < n; ++j) {
      if (j != i) {
        if (dist_mat.get_val(i, j) < min_val) min_val = dist_mat.get_val(i, j);
      }
    }
    mntd += min_val;
  }
  mntd *= 1.0 / n;
  return(mntd);
}


double calc_mpd_stat(const std::vector< std::array< size_t, 2 >>& edge,
                         const std::vector<double>& el) {
  auto dist_mat = dist_nodes_tri(edge, el);
  //int n = (el.size() + 2) * 0.5;
  //int max_pos = n * (n - 1) * 0.5;
  int max_pos = 0.125 * (el.size() * el.size()) + 0.25 * el.size();
  double mpd = std::accumulate(dist_mat.data_.begin(),
                               dist_mat.data_.begin() + max_pos, 0.0);
  mpd *= 1.0 / max_pos;
  return mpd;
}

double calc_psv_stat(const std::vector< std::array< size_t, 2 >>& edge,
                     const std::vector<double>& el) {

  auto dist_mat = dist_nodes_tri(edge, el);
  int n = (el.size() + 2) * 0.5;
  int max_pos = 0.125 * (el.size() * el.size()) + 0.25 * el.size();
  auto psv = 0.5 * std::accumulate(dist_mat.data_.begin(),
                              dist_mat.data_.begin() + max_pos, 0.0);

  psv *= 1.0 / (n * (n - 1));
  psv *= 2.0; // post hoc correction to match picante::psv output, because we
  // wrongly measure distance to most common ancestor per node.
  return(psv);
}

double calc_var_mpd_stat(const std::vector< std::array< size_t, 2 >>& edge,
                         const std::vector<double>& el) {
  auto dist_mat = dist_nodes_tri(edge, el);

  int max_pos = 0.125 * (el.size() * el.size()) + 0.25 * el.size();

  double s = 0.0;
  double s2 = 0.0;

  for (auto it = dist_mat.data_.begin(); it != dist_mat.data_.begin() + max_pos; ++it) {
    s  += (*it);
    s2 += (*it) * (*it);
  }

  double inv_max_pos = 1.0 / max_pos;

  return (s2 - (s * s) * inv_max_pos) * inv_max_pos;
}

#endif
