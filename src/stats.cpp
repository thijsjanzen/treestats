#include <Rcpp.h>

#include "beta.h"
#include "nltt.h"
#include "sackin.h"
#include "gamma.h"
#include "phylodiv.h"
#include "pigot_rho.h"


// [[Rcpp::export]]
double calc_beta_cpp(Rcpp::NumericMatrix in_table,
                     double upper_lim) {
  if (in_table.ncol() != 4) {
    Rcpp::stop("ltable needs four columns");
  }
  ltable ltab;
  for (int i = 0; i < in_table.nrow(); ++i) {
    std::array< float, 4 > row;
    for (int j = 0; j < in_table.ncol(); ++j) {
      row[j] = in_table(i,j);
    }
    ltab.push_back(row);
  }

  //  force_output("collected into vector\n");

  try {
    double output = calc_beta(ltab, -2, upper_lim);
    return output;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// [[Rcpp::export]]
float calc_sackin_cpp(Rcpp::NumericMatrix in_table,
                    std::string normalization) {

  std::vector< std::vector< float >> ltab;
  for (int i = 0; i < in_table.nrow(); ++i) {
    std::vector< float > row;
    for (int j = 0; j < in_table.ncol(); ++j) {
      row.push_back(in_table(i,j));
    }
    ltab.push_back(row);
  }

  //  force_output("collected into vector\n");

  try {
    float output = static_cast<float>(calc_sackin(ltab));

    if (normalization == "yule") {
      output = correct_yule(ltab, output);
    }
    if (normalization == "pda") {
      output = correct_pda(ltab, output);
    }

    return output;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// [[Rcpp::export]]
float calc_nltt_cpp(const Rcpp::NumericVector& brts_one,
                    const Rcpp::NumericVector& brts_two) {

try {
  std::vector<float> v1(brts_one.begin(), brts_one.end());
  std::vector<float> v2(brts_two.begin(), brts_two.end());

  auto nltt = calc_nltt(v1, v2);
  return nltt;
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}

// [[Rcpp::export]]
float calc_gamma_cpp(const Rcpp::NumericVector& brts_in) {
try {
  std::vector<float> brts(brts_in.begin(), brts_in.end());
  gamma_stat gamma_stat_object(brts);

  return gamma_stat_object.calc_gamma_stat();
} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}

// [[Rcpp::export]]
double calc_phylodiv_cpp(const Rcpp::List& phy,
                         double t,
                         double crown_age,
                         double extinct_acc) {
try {

  Rcpp::NumericMatrix edge = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  std::vector<double> el(edge_length.begin(), edge_length.end());
  std::vector< std::array<int, 2>> edges(edge.nrow());
  for (size_t i = 0; i < edge.nrow(); i++) {
    std::array<int, 2> to_add = {static_cast<int>(edge(i, 0)),
                                 static_cast<int>(edge(i, 1))};
    edges[i] = to_add;
  }

  phylo phylo_tree(edges, el);

  return calculate_phylogenetic_diversity(phylo_tree, t, crown_age, extinct_acc);

} catch(std::exception &ex) {
  forward_exception_to_r(ex);
} catch(...) {
  ::Rf_error("c++ exception (unknown reason)");
}
return NA_REAL;
}


// [[Rcpp::export]]
double calc_rho_cpp(const Rcpp::List& phy,
                    double crown_age) {
  try {

    Rcpp::NumericMatrix edge = phy["edge"];
    Rcpp::NumericVector edge_length = phy["edge.length"];

    std::vector<double> el(edge_length.begin(), edge_length.end());
    std::vector< std::array<int, 2>> edges(edge.nrow());
    for (size_t i = 0; i < edge.nrow(); i++) {
      std::array<int, 2> to_add = {static_cast<int>(edge(i, 0)),
                                   static_cast<int>(edge(i, 1))};
      edges[i] = to_add;
    }

    phylo phylo_tree(edges, el);
    rho pigot_rho(phylo_tree, crown_age);
    return pigot_rho.calc_pigot_rho();

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

struct bl {
  float e1, e2, EL;
};

//' function to calculate branching times of a tree
//' @param phy phylo object
//' @export
// [[Rcpp::export]]
std::vector< float > branching_times(const Rcpp::List& phy) {

  size_t Nnode = phy["Nnode"];
  size_t n = Nnode + 1;

  std::vector<size_t> interns(Nnode);

  std::vector<float> edge_length = phy["edge.length"];
  Rcpp::NumericMatrix edge = phy["edge"];

  size_t cnt = 0;
  for (size_t i = 0; i < edge_length.size(); ++i) {
    if (edge(i, 1) > n) {
      interns[cnt] = i;
      cnt++;
    }
  }

  std::vector<float> xx(Nnode);

  for (const auto& i : interns) {
    xx[ edge(i, 1) - n - 1 ] = xx[edge(i, 0) - n - 1] + edge_length[i];
  }

  int N = edge_length.size() - 1;
  float depth = xx[edge(N, 0) - n - 1] +  edge_length[N];
  for (auto& i : xx) {
    i = depth - i;
  }
  return xx;
}
