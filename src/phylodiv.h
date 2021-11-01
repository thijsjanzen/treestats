#ifndef phylo_div_h
#define phylo_div_h

struct branch {

  branch(float bd, int pl, int lab, float ext, float branch_length) :
  start_date(bd),
  parent_label(pl),
  label(lab),
  end_date(ext),
  bl(branch_length)
  {}

  float start_date;
  int parent_label;
  int label;
  float end_date;
  float bl;
};

struct phylo {

  const std::vector< std::array<int, 2>> edge;
  const std::vector< float > edge_length;

  phylo(const std::vector< std::array<int, 2>> e,
        const std::vector<float> el) : edge(e), edge_length(el) {}
};


float get_start_date(const std::vector<branch>& branchset,
                     int parent_label) {

  for (const auto& i : branchset) {
    if (i.label == parent_label) {
      return(i.end_date);
    }
  }
  // throw std::runtime_error("can not find parent");
  return 1e10f; // max MAX val possible
}

bool has_no_daughters(const std::vector<branch>& bs,
                      size_t parent_label) {
  // find any branches that share this parent
  for (const auto& i : bs) {
    if (i.parent_label == parent_label) {
      return false;
    }
  }
  return true;
}

std::vector< branch > remove_from_branchset(std::vector<branch> bs,
                                            size_t label) {

  size_t index = 0;
  for (; index < bs.size(); ++index) {
    if (bs[index].label == label)
      break;
  }
  size_t parent_label = bs[index].parent_label;
  // remove focal branch
  bs[index] = bs.back();
  bs.pop_back();

  // now check if parent has daughters
  if (has_no_daughters(bs, parent_label)) {
    bs = remove_from_branchset(bs, parent_label);
  }
  return bs;
}


std::vector< branch > create_branch_set(const phylo& phy,
                                        float max_t) {

  std::vector< branch > branchset;
  size_t crown = phy.edge[0][0];
  std::vector < std::array< float, 2 >> tip_times;
  float crown_age = 0.f;
  for (size_t i = 0; i < phy.edge.size(); ++i) {
    size_t parent_label = phy.edge[i][0];

    float start_date = 0.f;
    if (parent_label != crown) {
      start_date = get_start_date(branchset, parent_label);
    }

    if (start_date > max_t)
      continue;

    size_t own_label   = phy.edge[i][1];
    float bl = phy.edge_length[i];
    float end_date = start_date + bl;

    if (end_date > max_t) {
      end_date = max_t;
      bl = end_date - start_date;
    }
    if (own_label < crown) {
      tip_times.push_back({static_cast<float>(own_label), end_date});
      if (end_date > crown_age) crown_age = end_date;
    }

    branch new_branch(start_date, parent_label, own_label, end_date, bl);
    branchset.push_back(new_branch);
  }

  for (int i = 0; i < tip_times.size(); ++i) {
    //if (tip_times[i][1] < crown_age) { // extinct tip!
    if (crown_age - tip_times[i][1] > 1e-4f) {
      branchset = remove_from_branchset(branchset, tip_times[i][0]);
    }
  }

  return(branchset);
}

float calculate_phylogenetic_diversity(const phylo& phy,
                                       float t) {
  std::vector< branch > branchset = create_branch_set(phy, t);

  float pd = 0.f;
  for (auto it = branchset.cbegin(); it != branchset.cend(); ++it) {
    pd += it->bl;
  }
  return pd;
}


#endif
