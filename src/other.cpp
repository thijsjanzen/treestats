/// INDEP


#include <vector>
#include <array>
#include <Rcpp.h>

#include "util.h"
#include "beta.h"
#include "phylo2L.h"
#include "L2newick.h"
#include "newick_to_edge.h"
#include "avgladder.h"

// [[Rcpp::export]]
double calc_beta_cpp(const Rcpp::List& phy,
                     double upper_lim,
                     std::string algorithm,
                     double abs_tol,
                     double rel_tol) {

  try {
    Rcpp::NumericMatrix edge = phy["edge"];
    edge_table local_edge(edge.nrow());
    for (size_t i = 0; i < edge.nrow(); ++i) {
      local_edge[i] = {static_cast<size_t>(edge(i, 0)),
                       static_cast<size_t>(edge(i, 1))};
    }

    return calc_beta(local_edge, -2.0, upper_lim, algorithm, abs_tol, rel_tol);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}


// [[Rcpp::export]]
double calc_beta_ltable_cpp(const Rcpp::NumericMatrix& ltable,
                         double upper_lim,
                         std::string algorithm,
                         double abs_tol,
                         double rel_tol) {

  try {
    std::vector< std::array< double, 4 >> ltab(ltable.nrow());
    for (size_t i = 0; i < ltable.nrow(); ++i) {
      std::array< double, 4> row_entry = {ltable(i, 0), ltable(i, 1),
                                          ltable(i, 2), ltable(i, 3)};
      ltab[i] = row_entry;
    }

    return calc_beta(ltab, -2.0, upper_lim, algorithm, abs_tol, rel_tol);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}


//' Function to generate an ltable from a phy object. This function is a C++
//' implementation of the function DDD::phylo2L. An L table summarises a
//' phylogeny in a table with four columns, being: 1) time at which a species
//' is born, 2) label of the parent of the species, where positive and negative
//' numbers indicate whether the species belongs to the left or right crown
//' lineage, 3) label of the daughter species itself (again positive or
//' negative depending on left or right crown lineage), and the last column 4)
//' indicates the time of extinction of a species, or -1 if the species is
//' extant.
//' @param phy phylo object
//' @export
//' @examples simulated_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
//' ltable <- phylo_to_l(simulated_tree)
//' reconstructed_tree <- DDD::L2phylo(ltable)
//' par(mfrow=c(1, 2))
//' # trees should be more or less similar, although labels may not match, and
//' # rotations might cause (initial) visual mismatches
//' plot(simulated_tree)
//' plot(reconstructed_tree)
// [[Rcpp::export]]
Rcpp::NumericMatrix phylo_to_l(const Rcpp::List& phy) {
  const size_t ncol = 4;
  std::vector< std::array< double, ncol> > ltab = phylo_to_l_cpp(phy);

  size_t nrow = ltab.size();
  Rcpp::NumericMatrix out(nrow, ncol);

  for (size_t i = 0; i < ltab.size(); ++i) {
    for (size_t j = 0; j < ncol; ++j) {
      out(i, j) = ltab[i][j];
    }
  }
  return out;
}


std::vector< std::vector< double >> dist_nodes(const Rcpp::List& phy) {
  Rcpp::NumericMatrix edge = phy["edge"];
  Rcpp::NumericVector el   = phy["edge.length"];

  int n = 1 + edge.nrow() / 2;
  int m = n - 1;

  // code below is from the Ape package
  std::vector< size_t > e1(edge.nrow());
  std::vector< size_t > e2(edge.nrow());

  for (size_t i = 0; i < edge.nrow(); ++i) {
    e1[i] = edge(i, 0) - 1;
    e2[i] = edge(i, 1) - 1;
  }

  int i, j, k, a, d, NM = n + m, ROOT;
  double x;
  size_t N = e1.size();
  std::vector< std::vector<double>> D(NM, std::vector<double>(NM, 0.0));

  ROOT = e1[0]; d = e2[0]; /* the 2 nodes of the 1st edge */
  D[ROOT][d] = D[d][ROOT] = -el[0]; /* the 1st edge gives the 1st distance */

  /* go down along the edge matrix
   starting at the 2nd edge: */
  for (i = 1; i < N; i++) {
    a = e1[i]; d = e2[i]; x = el[i]; /* get the i-th nodes and branch length */
  D[a][d] = D[d][a] = -x;
  /* then go up along the edge matrix from the i-th edge
   to visit the nodes already visited and update the distances: */
  for (j = i - 1; j >= 0; j--) {
    k = e2[j];
    if (k == a) continue;
    D[k][d] = D[d][k] = D[a][k] - x;
  }
  if (k != ROOT)
    D[ROOT][d] = D[d][ROOT] = D[ROOT][a] - x;
  }
  return D;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix prep_lapl_spec(const Rcpp::List& phy) {
  std::vector< std::vector< double >> lapl_mat = dist_nodes(phy);
  Rcpp::NumericMatrix res(lapl_mat.size(), lapl_mat[0].size());

  for (size_t i = 0; i < lapl_mat.size(); ++i) {
    for (size_t j = 0; j < lapl_mat[i].size(); ++j) {
      res(i, j) = lapl_mat[i][j];
    }
    res(i, i) = - std::accumulate(lapl_mat[i].begin(), lapl_mat[i].end(), 0.0);
  }

  return res;
}


// [[Rcpp::export]]
double calc_mpd_cpp(const Rcpp::List& phy) {
  auto dist_mat = dist_nodes(phy);

  std::vector< std::string> tip_labels = phy["tip.label"];
  int n = tip_labels.size();
  double mpd = 0.0;
  size_t cnt = 0;


  for (size_t i = 0; i <n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      mpd += std::abs(dist_mat[i][j]); cnt++;
    }
  }
  mpd *= 1.0 / cnt;
  return(mpd);
}


// [[Rcpp::export]]
double calc_psv_cpp(const Rcpp::List& phy) {
  auto dist_mat = dist_nodes(phy);

  std::vector< std::string> tip_labels = phy["tip.label"];
  int n = tip_labels.size();
  double psv = 0.0;

  for (size_t i = 0; i <n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      psv += 0.5 * std::abs(dist_mat[i][j]);
    }
  }
  psv *= 1.0 / (n * (n - 1));
  psv *= 2.0; // post hoc correction to match picante::psv output, because we
              // wrongly measure distance to most common ancestor per node.
  return(psv);
}



double calc_J(const Rcpp::List& phy) {
  auto dist_mat = dist_nodes(phy);

  std::vector< std::string> tip_labels = phy["tip.label"];
  int n = tip_labels.size();
  double mpd = 0.0;
  size_t cnt = 0;


  for (size_t i = 0; i <n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      mpd += std::abs(dist_mat[i][j]); cnt++;
    }
  }
  return(mpd * 1.0 / (n * n));
}



// [[Rcpp::export]]
double calc_mntd_cpp(const Rcpp::List& phy) {
  auto dist_mat = dist_nodes(phy);

  std::vector< std::string> tip_labels = phy["tip.label"];
  int n = tip_labels.size();

  double mntd = 0.0;
  for (size_t i = 0; i < n; ++i) {
    dist_mat[i][i] = 1e6;
    double min_val = std::abs(dist_mat[i][0]);
    for (size_t j = 1; j < n; ++j) {
      if (std::abs(dist_mat[i][j]) < min_val) min_val = std::abs(dist_mat[i][j]);
    }
    mntd += min_val;
  }
  mntd *= 1.0 / n;
  return(mntd);
}


// [[Rcpp::export]]
double calc_var_mpd_cpp(const Rcpp::List& phy) {
  auto dist_mat = dist_nodes(phy);
  std::vector< std::string> tip_labels = phy["tip.label"];
  int n = tip_labels.size();
  double mpd = 0.0;
  size_t cnt = 0;
  for (size_t i = 0; i <n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      dist_mat[i][j] = std::abs(dist_mat[i][j]);
      mpd += dist_mat[i][j]; cnt++;
    }
  }
  mpd *= 1.0 / cnt;

  double var_mpd = 0.0;
  for (size_t i = 0; i <n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      var_mpd += (dist_mat[i][j] - mpd) * (dist_mat[i][j] - mpd);
    }
  }
  var_mpd *= 1.0 / cnt;
  return var_mpd;
}


// [[Rcpp::export]]
double avgLadder_cpp(const std::vector<long>& tree_edge) {
  try {
  return calc_ladder(tree_edge);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// [[Rcpp::export]]
std::string l_to_newick(const Rcpp::NumericMatrix& ltable_R,
                        bool drop_extinct) {
  auto ltable_cpp = convert_to_ltable(ltable_R);
  auto newick_string = ltable_to_newick(ltable_cpp, drop_extinct);
  return newick_string;
}





// [[Rcpp::export]]
double avgLadder_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto newick_str = l_to_newick(ltable_R, false);
  auto edge_table = newick_to_edge(newick_str);
  return avgLadder_cpp(edge_table);
}



