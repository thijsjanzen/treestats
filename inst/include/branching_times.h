#ifndef BRANCHING_TIMES_H
#define BRANCHING_TIMES_H

#include <vector>
#include <array>

std::vector< double > branching_times(const std::vector< std::array< size_t, 2>>& edge,
                                      const std::vector<double>& edge_length,
                                      size_t Nnode,
                                      size_t n) {

  std::vector<double> xx(Nnode, 0.0);

 // for (size_t i = 0; i < edge.size(); ++i) {
//    std::cerr << edge[i][0] << " " << edge[i][1] << " " << edge_length[i] << "\n";
//  }


//  std::cerr << " starting interns loop\n";
 /* size_t i = 0;
  for (const auto& j : edge) {
    if (j[1] > n) {
      auto target_index = j[1] - Nnode - 2;
      auto source_index = j[0] - Nnode - 2;
      std::cerr << j[1] << " " << target_index << " " << source_index << "\n";
      xx [ target_index ] = xx[ source_index ] + edge_length[i];
    }
    i++;
  }*/

 for (size_t i = 0; i < edge_length.size(); ++i) {
   if (edge[i][1] > n) {
     auto target_index = edge[i][1] - n - 1; // e2[i] - n, -2 because -1 of R, and -1 of n (R->CPP conversion)
     auto source_index = edge[i][0] - n - 1; // e1[i] - n
   //  std::cerr << edge[i][1] << " " << target_index << " " << source_index << "\n";
     xx [ target_index ] = xx[ source_index ] + edge_length[i];
   }
 }


  auto edge_index = edge[edge_length.size() - 1][0] - n - 1;

  //std::cerr << edge_index << "\n";

  double depth = xx[edge_index] + edge_length.back();

 // std::cerr << depth << "\n";
 // for (auto i : xx) {
 ////   std::cerr << i << " ";
 // } std::cerr << "\n";


  for (auto& i : xx) {
    i = depth - i;
  }
  return xx;
}

#endif
