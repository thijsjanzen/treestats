// Copyright 2022 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
#pragma once

#include <vector>
#include <array>


double calc_ladder(const std::vector< int >& tree_edge,
                   bool return_max) {
  struct node_entry {
    std::array< int, 2> daughters;
    size_t daughter_index = 0;

    void add_daughter(int index) {
      daughters[daughter_index] = index;
      daughter_index++;
    }
  };

  int max_node_val = *std::max_element(tree_edge.begin(), tree_edge.end());
  int root_no = tree_edge[0];
  for (size_t i = 2; i < tree_edge.size(); i+=2) {
    if (tree_edge[i] < root_no) root_no = tree_edge[i];
  }

  std::vector< node_entry > edge_mat(max_node_val + 1 - root_no);
  std::vector< int > tips(edge_mat.size(), 0);

  for (size_t i = 0; i < tree_edge.size(); i += 2) {
    int node_lab = tree_edge[i] - root_no;

    int other_lab = tree_edge[i + 1] - root_no;

    if (node_lab > static_cast<int>(edge_mat.size()) || node_lab < 0) {
      throw std::out_of_range("node_lab > edge_mat.size()");
    }

    edge_mat[node_lab].add_daughter(other_lab);
    if (other_lab < 0) {
      tips[node_lab]++;
    }
  }

  for (auto& i : tips) {
    if (i != 1) i = 0;
  }


  double store_val = 0.0;
  int count_val = 0;

  for (size_t i = 0; i < edge_mat.size(); ++i) {
    auto daughter1 = edge_mat[i].daughters[0];
    auto daughter2 = edge_mat[i].daughters[1];

    if (daughter1 > 0 && daughter1 > static_cast<int>(tips.size())) {
      throw std::out_of_range("daughter1 > tips.size()");
    }
    if (daughter2 > 0 && daughter2 > static_cast<int>(tips.size())) {
      throw std::out_of_range("daughter2 > tips.size()");
    }

    if ((daughter1 >= 0) && (tips[daughter1] == 1)) {
      tips[daughter1] += tips[i];
      tips[i] = 0;
    } else {
      if ((daughter2 >= 0)  && (tips[daughter2] == 1)) {
        tips[daughter2] += tips[i];
        tips[i] = 0;
      }
    }

    if (i > tips.size() || i < 0) {
      throw std::out_of_range("i > tips.size()");
    }

    if (tips[i] > 1)  {
      if (return_max == true) {
        if (tips[i] > store_val) store_val = tips[i];
      } else {
        store_val += tips[i];
      }
      count_val++;
    }
  }

  if (count_val > 0 && return_max == false) store_val *= 1.0 / count_val;

  return store_val;  // * 1.0 / count_val;
}
