#ifndef phylo_to_l_h
#define phylo_to_l_h

#include <Rcpp.h>
#include <vector>
#include <array>
#include <algorithm>

#include <thread>
#include <chrono>

#include <iostream>

#include "branching_times.h"

void force_output() {
  // std::this_thread::sleep_for(std::chrono::nanoseconds(100));
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}


size_t get_min_index(const std::vector< std::array<double, 6>>& localtab,
                     size_t col_index) {

  size_t index = 0;
  double min_val = localtab[index][col_index];
  for (size_t i = 1; i < localtab.size(); ++i) {
    if (localtab[i][col_index] < min_val) {
      min_val = localtab[i][col_index];
      index = i;
    }
  }
  return index;
}

bool parent_in_nodesindex(const std::vector< size_t >& nodesindex,
                          size_t parent) {
  bool is_present = false;
  for (const auto& i : nodesindex) {
    if (i == parent) {
      is_present = true;
      break;
    }
  }
  return is_present;
}


void remove_from_L(std::vector< std::array<double, 6>>& L,
                   size_t j) {
  std::swap(L[j], L.back());
  L.pop_back();
}


std::vector< std::array<double, 6>> get_realL(const std::vector< size_t >& nodesindex,
                                              std::vector< std::array<double, 6>> L) {

  std::vector< std::array<double, 6>> realL;

  while(true) {
    size_t j = get_min_index(L, 2);
    size_t daughter = L[j][2];
    size_t parent   = L[j][1];

    if (parent_in_nodesindex(nodesindex, parent)) {

      for (size_t index = 0; index < L.size(); ++index) {
        if (L[index][1] == parent) {
          L[index][1] = daughter;
        }
      }

      bool match_found = false;
      for (auto& i : L) {
        if (i[2] == parent) {
          i[5] = L[j][5];
          i[2] = daughter;
          remove_from_L(L, j);
          match_found = true;
          break;
        }
      }
      if (!match_found) {
        realL.push_back(L[j]);
        remove_from_L(L, j);
      }

    } else {
      realL.push_back(L[j]);
      remove_from_L(L, j);
    }

    if (L.empty()) {
      break;
    }
  }

  std::sort(realL.begin(), realL.end(), [&](const std::array< double, 6>& v1,
                        const std::array< double, 6>& v2) {
    return(v1[0] > v2[0]); // sort decreasing
  });

  return realL;
}


std::vector< std::array< double, 4> > phylo_to_l_cpp(const Rcpp::List& phy) {
//  Rcpp::Rcout << "getting brts\n"; force_output();
  std::vector< double > brts = branching_times(phy);

//  Rcpp::Rcout << "starting min_brts\n"; force_output();
  auto min_brts = *std::min_element(brts.begin(), brts.end());
  if (min_brts < 0) {
    for (auto& i : brts) {
      i += fabs(min_brts);
    }
  }

 // Rcpp::Rcout << "extracting phy properties\n"; force_output();
  size_t num_species = static_cast<size_t>(phy["Nnode"]) + 1;

  Rcpp::StringVector tiplabel = phy["tip.label"];
  Rcpp::NumericMatrix edge    = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  size_t num_tips = tiplabel.size();

  std::vector< long double > brt_preL(edge.nrow());
  long double min_brt_preL = 1e10;

 // Rcpp::Rcout << "starting brt_preL\n"; force_output();

  for (size_t i = 0; i < edge.nrow(); ++i) {
    auto index = edge(i, 0) - num_tips - 1; // -1 because 0 indexing
    brt_preL[i] = brts[index];
    if (brt_preL[i] < min_brt_preL) {
      min_brt_preL = brt_preL[i];
    }
  }

  if (min_brt_preL == 0.0) {
    long double correction = 0.0;
    for (size_t i = 0; i < edge_length.size(); ++i) {
      if (brt_preL[i] == 0.0) {
        if (edge_length[i] > correction) {
          correction = edge_length[i];
        }
      }
    }

    for (auto& i : brt_preL) {
      i += correction;
    }
  }

  std::vector< std::array<double, 6>> pre_Ltable(brt_preL.size());

  for (size_t i = 0; i < brt_preL.size(); ++i) {
    pre_Ltable[i][0] = brt_preL[i];
    pre_Ltable[i][1] = edge(i, 0);
    pre_Ltable[i][2] = edge(i, 1);
    pre_Ltable[i][3] = edge_length[i];
    pre_Ltable[i][4] = brt_preL[i] - edge_length[i];
  }

/*  std::cerr << "pre_Ltable: \n";
  for (auto i : pre_Ltable) {
    for (auto j : i) {
      std::cerr << j << " ";
    }
    std::cerr << '\n';
  }*/



  // pre_Ltable confirmed correct
  //
  //
//  Rcpp::Rcout << "starting eeindicator\n"; force_output();
  std::vector<double> eeindicator(edge_length.size(), 0);

  std::vector< size_t > extant_species_index;
  for (size_t i = 0; i < pre_Ltable.size(); ++i) {

    if (pre_Ltable[i][4] <= 1e-6) {
      extant_species_index.push_back(pre_Ltable[i][2]);

      size_t index = 0;
      for (; index < pre_Ltable.size(); ++index) {
        if (pre_Ltable[index][2] == pre_Ltable[i][2]) {
          break;
        }
      }
      if (index != pre_Ltable.size()) {
        eeindicator[index] = -1;
      }

    }
  }
 /* std::cerr << "extant_species_index: ";
  for (auto i : extant_species_index) {
    std::cerr << i << " ";
  } std::cerr << "\n";*/

  /*
   std::vector< size_t > tipsindex;
   for (size_t i = 1; i <= num_species; ++i) {
   tipsindex.push_back(i);
   }

   */

//  Rcpp::Rcout << "starting extant_species_index\n"; force_output();
  for (size_t i = 1; i <= num_species; ++i) {
    bool found = false;
    for (const auto& j : extant_species_index) {
      if (i == j) {
        found = true;
        break;
      }
    }
    if (!found) {
      auto extinct_index3_local = i;
      size_t index = 0;
      for (; index < pre_Ltable.size(); ++index) {
        if (pre_Ltable[index][2] == extinct_index3_local) {
          break;
        }
      }
      if (index != pre_Ltable.size()) {
        eeindicator[index] = pre_Ltable[index][4];
      }
    }
  }

  for (const auto& i : extant_species_index) {
    size_t index = 0;
    for (; index < pre_Ltable.size(); ++index) {
      if (pre_Ltable[index][2] == i) {
        break;
      }
    }
    if (index != pre_Ltable.size()) {
      eeindicator[index] = -1;
    }
  }

  // verified correct so far


 // Rcpp::Rcout << "starting eeindicator2\n"; force_output();

  for (size_t i = 0; i < eeindicator.size(); ++i) {
    pre_Ltable[i][5] = eeindicator[i];
  }


  std::sort(pre_Ltable.begin(), pre_Ltable.end(), [&](const std::array< double, 6>& v1,
                             const std::array< double, 6>& v2) {
    return(v1[0] > v2[0]); // sort decreasing
  });

//  Rcpp::Rcout << "starting nodesindex\n"; force_output();
  std::vector< size_t > nodesindex(edge.nrow());
  for (size_t i = 0; i < edge.nrow(); ++i) {
    nodesindex[i] = static_cast<size_t>(edge(i, 0));
  }

//  Rcpp::Rcout << "starting get_realL\n"; force_output();
  std::vector< std::array<double, 6>> realL = get_realL(nodesindex,
                                                        pre_Ltable);

  // verified correct so far
  //
  //
 // Rcpp::Rcout << "starting creation L\n"; force_output();
  std::vector< std::array< double, 4> > L( realL.size() );

  for (size_t i = 0; i < realL.size(); ++i) {
    L[i][0] = realL[i][0];
    L[i][1] = realL[i][1];
    L[i][2] = i + 1;
    L[i][3] = realL[i][5];

    size_t index = 0;
    for (; index < realL.size(); ++index) {
      if (realL[index][2] == realL[i][1]) {
        break;
      }
    }
    if (index != realL.size()) {
      L[i][1] = index + 1;
    }
  }

  L[0][1] = 0;
  L[0][2] = -1;
  L[1][1] = -1;

  for (size_t i = 1; i < L.size(); ++i) {
    if (L[i - 1][2] < 0) {
      auto ref = fabs(L[i - 1][2]);
      for (auto & j : L) {
        if (j[1] == ref) {
          j[1] = L[i - 1][2];
          j[2] *= -1;
        }
      }
    }
  }
//  Rcpp::Rcout << "done L\n"; force_output();
  return L;
}




#endif
