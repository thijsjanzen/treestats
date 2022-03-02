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


#endif
