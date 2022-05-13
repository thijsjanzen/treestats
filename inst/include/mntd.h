#ifndef mntd_h
#define mntd_h

#include <vector>
#include <array>
#include "dist_nodes.h"

using ltable = std::vector< std::array<double, 4>>;

double calc_mntd_ltable(const ltable& ltable_) {
  std::vector<double> dist(ltable_.size() + 1, -1);

  for (const auto& i : ltable_) {
    auto parent = std::abs(i[1]);
    auto daughter = std::abs(i[2]);
    auto dist_to_nearest_taxon = 2 * i[0];
    if (i[3] != -1) {
      dist_to_nearest_taxon = i[0] + (i[0] - i[3]);
    }
    dist[daughter] = dist_to_nearest_taxon;
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
    double min_val = std::abs(dist_mat[i][0]);
    for (size_t j = 1; j < n; ++j) {
      if (std::abs(dist_mat[i][j]) < min_val) min_val = std::abs(dist_mat[i][j]);
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

  for (size_t i = 0; i <n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      mpd += std::abs(dist_mat[i][j]); cnt++;
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
      psv += 0.5 * std::abs(dist_mat[i][j]);
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
      dist_mat[i][j] = std::abs(dist_mat[i][j]);
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
