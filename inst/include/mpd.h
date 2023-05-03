namespace mpd_tree {


struct node {
  node* daughter1 = nullptr;
  node* daughter2 = nullptr;
  size_t L;
  size_t R;
  double bl_R;
  double bl_L;

  node() {
    daughter1 = nullptr;
    daughter2 = nullptr;
    L = 0;
    R = 0;
    bl_R = -1.0;
    bl_L = -1.0;
  }

  size_t update_num_tips() {
    if (daughter1 && !daughter2) {
      L = daughter1->update_num_tips();
    }
    if (daughter1 && daughter2) {
      L =  daughter1->update_num_tips();
      R =  daughter2->update_num_tips();
    }

    return L + R;
  }
};

class phylo_tree {
public:
  explicit phylo_tree(const std::vector< int >& tree_edge,
                      const std::vector<double>& edge_length) {
    int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
    tree_size = root_no - 1;

    tree.resize(tree_edge.size() / 2 - root_no + 2);

    for (size_t i = 0; i < tree_edge.size(); i += 2) {
      int index    = static_cast<int>(tree_edge[i]) - root_no;
      int d1_index = static_cast<int>(tree_edge[i + 1]) - root_no;

      if (d1_index < 0) {
        tree[index].L == 0 ? tree[index].L = 1 : tree[index].R = 1;
      } else {
        !tree[index].daughter1 ?
         tree[index].daughter1 = &tree[d1_index] :
         tree[index].daughter2 = &tree[d1_index];
      }

      int el_index = i / 2;

      tree[index].bl_L < 0.0 ? tree[index].bl_L = edge_length[el_index] :
                               tree[index].bl_R = edge_length[el_index];
    }

    tree[0].update_num_tips();
  }

  double calculate_mpd() {
    int N = tree_size;
    double mpd = 0.0;
    for (const auto& i : tree) {
      int l = i.L;
      int r = i.R;
      auto L_bl = i.bl_L;
      auto R_bl = i.bl_R;

      std::cerr << l << " " << L_bl << " " << r << " " << R_bl << "\n";

      mpd += L_bl * (l * (N - l));
      mpd += R_bl * (r * (N - r));
    }

    // std::cerr << N << " " << mpd << "\n";

    mpd *= 2.0 / (N * (N - 1));
    return(mpd);
  }

private:
  std::vector< node > tree;
  int tree_size = 0;
};


} // namespace mpd_tree



/*
 struct mpd_tree {
 mpd_tree(const std::vector< int >& tree_edge,
 const std::vector<double>& edge_length) : branch_length(edge_length){
 // we need to construct a bl, descendants matrix
 int root_no = 2 + static_cast<int>(0.25 * tree_edge.size());
 num_descendants = std::vector<int>(branch_length.size());
 for (int i = 1; i < tree_edge.size(); ++i) {
 if (tree_edge[i] < root_no) {
 num_descendants[i / 2] = 1;
 }
 }
 for (int i = 0; i < tree_edge.size(); i += 2) {
 int parent = tree_edge[i];
 int offspring = tree_edge[i + 1];
 num_descendants[parent] += num_descendants[offspring];
 }

 }

 private:
 const std::vector<double> branch_length;
 std::vector<int> num_descendants;
 };
 */
