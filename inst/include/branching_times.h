#ifndef BRANCHING_TIMES_H
#define BRANCHING_TIMES_H

#include <vector>
#include <array>

std::vector< double > branching_times(const std::vector< std::array< size_t, 2>>& edge,
                                      const std::vector<double>& edge_length,
                                      size_t Nnode) {

  std::vector<double> xx(Nnode, 0.0);

  size_t i = 0;
  for (const auto& j : edge) {
    if (j[1] > Nnode + 1) {
      auto target_index = j[1] - Nnode - 2;
      auto source_index = j[0] - Nnode - 2;
      xx [ target_index ] = xx[ source_index ] + edge_length[i];
    }
    i++;
  }

  auto edge_index = edge[edge_length.size() - 1][0] - Nnode - 2;

  double depth = xx[edge_index] + edge_length.back();
  for (auto& i : xx) {
    i = depth - i;
  }
  return xx;
}

#endif
