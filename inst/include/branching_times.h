#ifndef BRANCHING_TIMES_H
#define BRANCHING_TIMES_H

#include <vector>
#include <array>

std::vector< double > branching_times(const std::vector< std::array< size_t, 2>>& edge,
                                      const std::vector<double> edge_length,
                                      size_t Nnode) {

  size_t n = Nnode + 1;

  std::vector<double> xx(Nnode, 0.0);

  for (size_t i = 0; i < edge_length.size(); ++i) {
    if (edge[i][1] > n) {
      auto target_index = edge[i][1] - n - 1; // -1 because of R to C++ indexing
      auto source_index = edge[i][0] - n - 1;

      xx[ target_index  ] = xx[ source_index] + edge_length[i];
    }
  }

  auto edge_index = edge[edge_length.size() - 1][0] - n - 1;

  double depth = xx[edge_index] +  edge_length.back();
  for (auto& i : xx) {
    i = depth - i;
  }
  return xx;
}

#endif
