#ifndef pigot_rho_h
#define pigot_rho_h

#include "phylodiv.h"

double calc_rho(const std::vector<double>& brts) {
  size_t n1 = 2;
  size_t n3 = 1 + brts.size();
  double mid_point = brts[0] / 2;

  auto mid_point_entry = std::lower_bound(brts.begin(), brts.end(), mid_point,
                                          std::greater<double>());
  size_t index_mid_point = std::distance(brts.begin(), mid_point_entry);

  size_t n2 = 1 + index_mid_point;

  double  r1 = (log(n2) - log(n1)) / mid_point;
  double  r2 = (log(n3) - log(n2)) / mid_point;

  return ((r2 - r1) / (r1 + r2));
}


struct rho {

  rho (const phylo& phy, double crown_age) {
    branchset = create_branch_set(phy, 1e12, crown_age, 1e-6);
    crown_age_ = crown_age;
    mid_point_ = crown_age_ * 0.5;
    extant_lineages = count_extant_lineages(phy);
  }


  size_t num_lineages(double t) {
    size_t num_lin = 0;
    for (const auto& i : branchset) {
      if (i.start_date < t && i.end_date >= t) num_lin++;
    }
    return num_lin;
  }

  double calc_pigot_rho() {

    size_t n1 = 2;
    size_t n2 = num_lineages(mid_point_);
    size_t n3 = extant_lineages;

    double r1 = (log(n2) - log(n1)) / mid_point_;
    double r2 = (log(n3) - log(n2)) / mid_point_;
    return (r2 - r1) / (r1 + r2);
  }

  size_t count_extant_lineages(const phylo& phy) {
    double crown = phy.edge[0][0];
    size_t extant_lin = 0;
    for (const auto& i : phy.edge) {
      if (i[1] < crown) extant_lin++;
    }
    return extant_lin;
  }

private:
  std::vector< branch > branchset;
  double crown_age_;
  double mid_point_;
  size_t extant_lineages;
};





#endif
