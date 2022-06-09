#ifndef mntd_h
#define mntd_h

#include <vector>
#include <array>
#include "dist_nodes.h"

using ltable = std::vector< std::array<double, 4>>;


template <typename T>
struct lower_tri {
  lower_tri(size_t n) : n_(n) {
    data_.resize( (n_ * (n_ - 1)) * 0.5);
  }

  void set_val(int i, int j, T val) {
    if (i == j) return; // do nothing, these values don't exist.

    if (i < j) data_[i * n_ + j] = val;
    if (i > j) data_[j * n_ + i] = val;
  }

  /*T get_val(int i, int j) {
    if (i == j) return static_cast<T>(0.0);

    if (i < j) return data_[i * n_ + j];

    return data_[j * n_ + i];
  }*/

  T get_val(int i, int j) {

   /* T out = 0.0;

    if (i < j) out = data_[i * n_ + j];
    if (j > i) out = data_[j * n_ + i];

    std::cerr << i << " " << j << " " << out << "\n";*/

    if (i == j) {
      std::cerr << i << " " << j << " " << -5.0 << "\n";
    }

    if (i < j) {
      std::cerr << i << " " << j << " " << i * n_ + j << "\n";
      return data_[i * n_ + j];
    }

    std::cerr << i << " " << j << " " << j * n_ + i << "\n";

    return data_[j * n_ + i];
  }

  std::vector<T> data_;
  size_t n_;
};

lower_tri<double> dist_nodes_tri(const std::vector< std::array< size_t, 2 >>& edge,
                         const std::vector<double>& el) {
  int n = 1 + edge.size() / 2;
  int m = n - 1;
  auto nm = n + m;
  static double max_s = 46340; // floor(sqrt(2^31 - 1))
  if (nm > max_s) {
    // std::cerr << n << " " << m << " " << nm << " " << max_s << "\n";
    throw std::runtime_error("tree too big");
  }
  // code below is from the Ape package
  std::vector< size_t > e1(edge.size());
  std::vector< size_t > e2(edge.size());

  for (size_t i = 0; i < edge.size(); ++i) {
    e1[i] = edge[i][0] - 1;
    e2[i] = edge[i][1] - 1;
  }

  int i, j, k, a, d, NM = n + m, ROOT;
  double x;
  size_t N = e1.size();
  lower_tri<double> D(NM);

  ROOT = e1[0]; d = e2[0]; /* the 2 nodes of the 1st edge */

  D.set_val(d, ROOT, el[0]);
 // D.set_val(ROOT, d, -el[0]);

  /* go down along the edge matrix
   starting at the 2nd edge: */
  for (i = 1; i < N; i++) {
      a = e1[i]; d = e2[i]; x = el[i]; /* get the i-th nodes and branch length */
    //D(a, d) = D(d, a) = -x;
      D.set_val(a, d, x);
    /* then go up along the edge matrix from the i-th edge
     to visit the nodes already visited and update the distances: */
    for (j = i - 1; j >= 0; j--) {
      k = e2[j];
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
  auto dist_mat = dist_nodes(edge, el);

  int n = (el.size() + 2) * 0.5;

  double mntd = 0.0;
  for (size_t i = 0; i < n; ++i) {
    dist_mat[i][i] = 1e6;
    double min_val = dist_mat[i][0];
    for (size_t j = 1; j < n; ++j) {
      if (dist_mat[i][j] < min_val) min_val = dist_mat[i][j];
    }
    mntd += min_val;
  }
  mntd *= 1.0 / n;
  return(mntd);
}

double calc_mpd_stat(const std::vector< std::array< size_t, 2 >>& edge,
                     const std::vector<double>& el) {
  int n = (el.size() + 2) * 0.5;
  auto dist_mat = dist_nodes(edge, el);

  double mpd = 0.0;
  size_t cnt = 0;

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      mpd += dist_mat[i][j]; cnt++;
    }
  }
  mpd *= 1.0 / cnt;
  return mpd;
}

double calc_mpd_stat_tri(const std::vector< std::array< size_t, 2 >>& edge,
                         const std::vector<double>& el) {
  int n = (el.size() + 2) * 0.5;
  auto dist_mat = dist_nodes_tri(edge, el);
  double mpd = 0.0;
  size_t cnt = 0;
  std::cerr << "starting collection\n";
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      mpd += dist_mat.get_val(i, j); cnt++;
    }
  }

  mpd *= 1.0 / cnt;
  return mpd;
}

double calc_psv_stat(const std::vector< std::array< size_t, 2 >>& edge,
                     const std::vector<double>& el) {

  auto dist_mat = dist_nodes(edge, el);
  int n = (el.size() + 2) * 0.5;

  double psv = 0.0;


  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      psv += 0.5 * dist_mat[i][j];
    }
  }
  psv *= 1.0 / (n * (n - 1));
  psv *= 2.0; // post hoc correction to match picante::psv output, because we
  // wrongly measure distance to most common ancestor per node.
  return(psv);
}

double calc_var_mpd_stat(const std::vector< std::array< size_t, 2 >>& edge,
                         const std::vector<double>& el) {
  auto dist_mat = dist_nodes(edge, el);
  int n = (el.size() + 2) * 0.5;

  double mpd = 0.0;
  size_t cnt = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      mpd += dist_mat[i][j]; cnt++;
    }
  }
  mpd *= 1.0 / cnt;

  double var_mpd = 0.0;
  for (size_t i = 0; i <n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      var_mpd += (dist_mat[i][j] - mpd) * (dist_mat[i][j] - mpd);
    }
  }
  var_mpd *= 1.0 / cnt;
  return var_mpd;
}

#endif
