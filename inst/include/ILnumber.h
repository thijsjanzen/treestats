#ifndef IL_H
#define IL_H

size_t calc_IL(const std::vector< long >& tree_edge) {

  int root_no = static_cast<int>(tree_edge.front());
  std::vector< int > nodes(tree_edge.size() / 2, 0);

  for (size_t i = 0; i < tree_edge.size(); i += 2 ) {
    if (tree_edge[i + 1] < root_no) {
      int index    = static_cast<int>(tree_edge[i]) - root_no;
      nodes[index]++;
    }
  }

  size_t num_IL = 0;
  for (const auto& i : nodes) {
    if (i == 1) num_IL++;
  }
  return(num_IL);
}

double calc_ladder2(const std::vector< long >& tree_edge) {

  auto max_node_val = *std::max_element(tree_edge.begin(), tree_edge.end());
  auto root_no = tree_edge[0];
  std::vector< std::array< size_t, 3 >> edge_mat(max_node_val + 1, {0, 0, 0});
  for (size_t i = 0; i < tree_edge.size(); i += 2 ) {
    auto node_lab = tree_edge[i];
    if (edge_mat[node_lab][0] == 0) {
      edge_mat[node_lab][0] = tree_edge[i+1];
    } else {
      edge_mat[node_lab][1] = tree_edge[i+1];
    }

    if (tree_edge[i + 1] < root_no) {
      edge_mat[node_lab][2]++;
    }
  }
  for (auto& i : edge_mat) {
    if (i[2] != 1) i[2] = 0;
  }

  for (auto& i : edge_mat) {
    auto daughter1 = i[0];
    auto daughter2 = i[1];
    if (edge_mat[daughter1][2] == 1) {
        edge_mat[daughter1][2] += i[2];
        i[2] = 0;
    } else {
      if (edge_mat[daughter2][2] == 1) {
        edge_mat[daughter2][2] += i[2];
        i[2] = 0;
      }
    }
  }
  double mean_val = 0.0;
  int count_val = 0;
  for (const auto& i : edge_mat) {
    if (i[2] > 1) {
      mean_val += i[2];
      count_val++;
    }
  }
  return mean_val * 1.0 / count_val;
}

double calc_ladder(const std::vector< long >& tree_edge) {
  struct node_entry {
    std::array< size_t, 2> daughters = {0, 0};
    size_t daughter_index = 0;
    size_t num_tips = 0;

    void add_daughter(size_t index) {
      daughters[daughter_index] = index;
      daughter_index++;
    }
  };

  auto max_node_val = *std::max_element(tree_edge.begin(), tree_edge.end());
  auto root_no = tree_edge[0];
  std::vector< node_entry > edge_mat(max_node_val + 1 - root_no);
  for (size_t i = 0; i < tree_edge.size(); i += 2 ) {
    auto node_lab = tree_edge[i] - root_no;

    auto other_lab = tree_edge[i+1] - root_no;
    edge_mat[node_lab].add_daughter(other_lab);
    if (other_lab < 0) edge_mat[node_lab].num_tips++;
  }

  for (auto& i : edge_mat) {
    if (i.num_tips != 1) i.num_tips = 0;
  }


  double mean_val = 0.0;
  int count_val = 0;

  for (auto& i : edge_mat) {

    auto daughter1 = i.daughters[0];
    auto daughter2 = i.daughters[1];
    if (edge_mat[daughter1].num_tips == 1) {
      edge_mat[daughter1].num_tips += i.num_tips;
      i.num_tips = 0;
    } else {
      if (edge_mat[daughter2].num_tips == 1) {
        edge_mat[daughter2].num_tips += i.num_tips;
        i.num_tips = 0;
      }
    }
    if (i.num_tips > 1)  {
      mean_val += i.num_tips;
      count_val++;
    }
  }

  return mean_val * 1.0 / count_val;
}


#endif
