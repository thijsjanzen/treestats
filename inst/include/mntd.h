#ifndef mntd_h
#define mntd_h

using ltable = std::vector< std::array<double, 4>>;

double calc_mntd_ltable(const ltable& ltable_) {
  std::vector<double> dist(ltable_.size() + 1, -1);

  for (const auto& i : ltable_) {
    auto parent = std::abs(i[1]);
    auto daughter = std::abs(i[2]);
    auto dist_to_nearest_taxon = 2 * std::abs(i[0]);
    if (i[3] != -1) {
      dist_to_nearest_taxon -= i[3]; // for extinct branches.
    }
    dist[daughter] = dist_to_nearest_taxon;
    if (dist[parent] > 0) {
      if (dist_to_nearest_taxon < dist[parent]) {
        dist[parent] = dist_to_nearest_taxon;
      }
    }
  }
  dist[0]= 0.0;
  auto sum_dist = std::accumulate(dist.begin(), dist.end(), 0.0);
  return sum_dist * 1.0 / ltable_.size();
}



#endif
