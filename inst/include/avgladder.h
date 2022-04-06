#ifndef AVGLADDER_H
#define AVGLADDER_H

#include <vector>
#include <array>


double calc_ladder(const std::vector< long >& tree_edge) {
  struct node_entry {
    std::array< int, 2> daughters;
    size_t daughter_index = 0;

    void add_daughter(int index) {
      daughters[daughter_index] = index;
      daughter_index++;
    }
  };

  long max_node_val = *std::max_element(tree_edge.begin(), tree_edge.end());
  int root_no = 2 + 0.25 * tree_edge.size();

  std::vector< node_entry > edge_mat(max_node_val + 1 - root_no);
  std::vector< int > tips(edge_mat.size(), 0);

  for (size_t i = 0; i < tree_edge.size(); i += 2 ) {
    int node_lab = tree_edge[i] - root_no;

    int other_lab = tree_edge[i + 1] - root_no;

    if (node_lab > edge_mat.size() || node_lab < 0) {
      for (size_t j = 0; j < tree_edge.size(); j +=2) {
        std::cerr << tree_edge[j] << " " << tree_edge[j + 1] << "\n";
      } std::cerr << "\n";
      std::cerr << node_lab << " " << other_lab << " " << root_no << "\n";

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


  double mean_val = 0.0;
  int count_val = 0;

  for (size_t i = 0; i < edge_mat.size(); ++i) {

    auto daughter1 = edge_mat[i].daughters[0];
    auto daughter2 = edge_mat[i].daughters[1];

    if (daughter1 > 0 && daughter1 >  tips.size()) {
      throw std::out_of_range("daughter1 > tips.size()");
    }
    if (daughter2 > 0 && daughter2 > tips.size()) {
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
      mean_val += tips[i];
      count_val++;
    }
  }

  if (count_val > 0) mean_val *= 1.0 / count_val;

  return mean_val; // * 1.0 / count_val;
}


#endif
