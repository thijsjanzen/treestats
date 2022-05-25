/// INDEP


#include <vector>
#include <array>
#include <Rcpp.h>

#include "dist_nodes.h"
#include "util.h"
#include "beta.h"
#include "phylo2L.h"
#include "L2newick.h"
#include "L2newick_ed.h"
#include "newick_to_edge.h"
#include "avgladder.h"
#include "mntd.h"

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




// [[Rcpp::export]]
double calc_mpd_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return calc_mpd_stat(edge, el);
}

// [[Rcpp::export]]
double calc_psv_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return calc_psv_stat(edge, el);
}

// [[Rcpp::export]]
double calc_J_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  auto mpd = calc_mpd_stat(edge, el);
  int n = (el.size() + 2) * 0.5;

  return mpd * 1.0 / n;
}

// [[Rcpp::export]]
double calc_mntd_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return calc_mntd_stat(edge, el);
}

// [[Rcpp::export]]
double calc_mntd_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
    auto ltab = convert_to_ltable(ltable_R);
    return calc_mntd_ltable(ltab);
}

// [[Rcpp::export]]
double calc_var_mpd_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return calc_var_mpd_stat(edge, el);
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
std::string l_to_newick_ed(const Rcpp::NumericMatrix& ltable_R, const double t,
                        bool drop_extinct) {
  auto ltable_cpp = convert_to_ltable(ltable_R);
  auto newick_string = ltable_to_newick_ed(ltable_cpp, t, drop_extinct);
  return newick_string;
}

// [[Rcpp::export]]
double avgLadder_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto newick_str = l_to_newick(ltable_R, false);
  auto edge_table = newick_to_edge(newick_str);
  return avgLadder_cpp(edge_table);
}
