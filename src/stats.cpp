#include <Rcpp.h>

#include <cstring>

#include "beta.h"
#include "nltt.h"
#include "sackin.h"
#include "gamma.h"
#include "phylodiv.h"
#include "pigot_rho.h"
#include "phylo2L.h"
#include "colless.h"
#include "L2newick.h"

#include "crown_age.h"

using ltable = std::vector< std::array<double, 4>>;

auto convert_to_ltable(const Rcpp::NumericMatrix& mat_in) {
  ltable out(mat_in.nrow());

  for (size_t i = 0; i < mat_in.nrow(); ++i) {
    std::array<double, 4> row_entry = {mat_in(i, 0), mat_in(i, 1),
                                       mat_in(i, 2), mat_in(i, 3) };
    out[i] = row_entry;
  }
  return out;
}

std::vector<double> branching_times_from_ltable(const Rcpp::NumericMatrix& mat_in) {
  std::vector<double> out(mat_in.nrow() - 1);
  for (size_t i = 1; i < mat_in.nrow(); ++i) {
    out[i - 1] = mat_in(i, 0);
  }
  return out;
}


// [[Rcpp::export]]
double calc_beta_cpp(const Rcpp::List& phy,
                     double upper_lim,
                     std::string algorithm,
                     double abs_tol,
                     double rel_tol) {

  try {
    Rcpp::NumericMatrix edge = phy["edge"];
    std::vector< std::array< size_t, 2 >> local_edge(edge.nrow());
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


// [[Rcpp::export]]
double calc_colless_cpp(const Rcpp::List phy,
                        std::string normalization) {

  Rcpp::NumericMatrix edge = phy["edge"];
  int num_tips = 1 + edge.nrow() / 2;
  int num_nodes = num_tips - 1;

  std::vector< int > parents = std::vector< int >(num_tips + num_nodes + 1, -1);

  for (size_t i = 0; i < edge.nrow(); ++i) {
    parents[ edge(i, 1) ] = edge(i, 0); // store parent ID for each tip / node
  }

  colless_stat s(parents, num_tips);

  double output = static_cast<double>(s.calc_colless());

  if (normalization == "yule") {
    output = s.correct_yule(output);
  }
  if (normalization == "pda") {
    output = s.correct_pda(output);
  }
  return output;

}

// [[Rcpp::export]]
double calc_colless_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                                std::string normalization) {

  auto l_in_cpp = convert_to_ltable(l_from_R);

  colless_stat_ltable s(l_in_cpp);

  double output = static_cast<double>(s.calc_colless());

  if (normalization == "yule") {
    output = s.correct_yule(output);
  }
  if (normalization == "pda") {
    output = s.correct_pda(output);
  }
  return output;
}

// [[Rcpp::export]]
double calc_blum_ltable_cpp(const Rcpp::NumericMatrix& ltab_in) {
  auto local_ltab = convert_to_ltable(ltab_in);
  sackin_stat_ltab s(local_ltab);

  return s.calc_blum();
}

// [[Rcpp::export]]
double calc_blum_cpp(const Rcpp::List phy) {

  Rcpp::NumericMatrix edge = phy["edge"];
  int num_tips = 1 + edge.nrow() / 2;
  int num_nodes = num_tips - 1;

  std::vector< int > parents = std::vector< int >(num_tips + num_nodes + 1, -1);

  for (size_t i = 0; i < edge.nrow(); ++i) {
    parents[ edge(i, 1) ] = edge(i, 0); // store parent ID for each tip / node
  }

  sackin_stat2 s(parents, num_tips);

  return s.calc_blum();
}

// [[Rcpp::export]]
double calc_sackin_cpp(const Rcpp::List phy,
                       const Rcpp::String& normalization) {

  Rcpp::NumericMatrix edge = phy["edge"];
  int num_tips = 1 + edge.nrow() / 2;
  int num_nodes = num_tips - 1;

  std::vector< int > parents = std::vector< int >(num_tips + num_nodes + 1, -1);

  for (size_t i = 0; i < edge.nrow(); ++i) {
    parents[ edge(i, 1) ] = edge(i, 0); // store parent ID for each tip / node
  }

  sackin_stat2 s(parents, num_tips);

  double output = static_cast<double>(s.calc_sackin());

  if (normalization == "yule") {
    Rcpp::StringVector tip_label = phy["tip.label"];
    size_t n = tip_label.size();
    output = s.correct_yule(n, output);
  }
  if (normalization == "pda") {
    Rcpp::StringVector tip_label = phy["tip.label"];
    size_t n = tip_label.size();
    output = s.correct_pda(n, output);
  }

  return output;
}


// [[Rcpp::export]]
double calc_sackin_ltable_cpp(const Rcpp::NumericMatrix& ltab,
                       const Rcpp::String& normalization) {

  auto local_ltab = convert_to_ltable(ltab);
  sackin_stat_ltab s(local_ltab);

  double output = static_cast<double>(s.calc_sackin());

  if (normalization == "yule") {
    output = s.correct_yule(output);
  }
  if (normalization == "pda") {
    output = s.correct_pda(output);
  }

  return output;
}


// [[Rcpp::export]]
double calc_nltt_ltable_cpp(const Rcpp::NumericMatrix& ltab1,
                            const Rcpp::NumericMatrix& ltab2) {

  auto brts_one = branching_times_from_ltable(ltab1);
  auto brts_two = branching_times_from_ltable(ltab2);
  std::sort(brts_one.begin(), brts_one.end(), std::greater<double>());
  std::sort(brts_two.begin(), brts_two.end(), std::greater<double>());
  for (auto& i : brts_one) {
    i *= -1;
  }
  for (auto& i : brts_two) {
    i *= -1;
  }
  brts_one.push_back(0.0);
  brts_two.push_back(0.0);
  auto nltt = calc_nltt(brts_one, brts_two);
  return nltt;
}

// [[Rcpp::export]]
double calc_nltt_cpp(const Rcpp::List& phy1,
                     const Rcpp::List& phy2) {

  std::vector<double> brts_one = branching_times(phy1);
  std::vector<double> brts_two = branching_times(phy2);
  std::sort(brts_one.begin(), brts_one.end(), std::greater<double>());
  std::sort(brts_two.begin(), brts_two.end(), std::greater<double>());
  for (auto& i : brts_one) {
    i *= -1;
  }
  for (auto& i : brts_two) {
    i *= -1;
  }
  brts_one.push_back(0.0);
  brts_two.push_back(0.0);
  auto nltt = calc_nltt(brts_one, brts_two);
  return nltt;
}


// [[Rcpp::export]]
double calc_gamma_cpp(const Rcpp::List& phy) {
  std::vector<double> brts = branching_times(phy);
  return calc_gamma(brts);
}

// [[Rcpp::export]]
double calc_gamma_ltable_cpp(const Rcpp::NumericMatrix& ltab_in) {
  std::vector<double> brts = branching_times_from_ltable(ltab_in);
  return calc_gamma(brts);
}

// [[Rcpp::export]]
double calc_phylodiv_cpp(const Rcpp::List& phy,
                         double t,
                         double extinct_acc) {
  try {

    Rcpp::NumericMatrix edge = phy["edge"];
    Rcpp::NumericVector edge_length = phy["edge.length"];

    std::vector<double> el(edge_length.begin(), edge_length.end());
    std::vector< std::array<size_t, 2>> edges(edge.nrow());
    for (size_t i = 0; i < edge.nrow(); i++) {
      std::array<size_t, 2> to_add = {static_cast<size_t>(edge(i, 0)),
                                   static_cast<size_t>(edge(i, 1))};
      edges[i] = to_add;
    }

    double crown_age = calc_crown_age(edges, el); // ignore root edge
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
double calc_rho_complete_cpp(const Rcpp::List& phy) {
  Rcpp::NumericMatrix edge = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  std::vector<double> el(edge_length.begin(), edge_length.end());
  std::vector< std::array<size_t, 2>> edges(edge.nrow());
  for (size_t i = 0; i < edge.nrow(); i++) {
    std::array<size_t, 2> to_add = {static_cast<size_t>(edge(i, 0)),
                                    static_cast<size_t>(edge(i, 1))};
    edges[i] = to_add;
  }

  double crown_age = calc_crown_age(edges, el);
  phylo phylo_tree(edges, el);
  rho pigot_rho(phylo_tree, crown_age);
  return pigot_rho.calc_pigot_rho();
}

// [[Rcpp::export]]
double calc_rho_cpp(const Rcpp::List& phy) {

  size_t num_nodes = static_cast<size_t>(phy["Nnode"]);

  if (num_nodes < 200) {
    return calc_rho_complete_cpp(phy);
  }

  auto brts = branching_times(phy);
  return calc_rho(brts);
}

// [[Rcpp::export]]
double calc_rho_ltable_cpp(const Rcpp::NumericMatrix& ltab) {

  auto brts = branching_times_from_ltable(ltab);
  return calc_rho(brts);
}


// [[Rcpp::export]]
double calc_crown_age_cpp(const Rcpp::List& phy) {
  Rcpp::NumericMatrix edge = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  std::vector<double> el(edge_length.begin(), edge_length.end());
  std::vector< std::array<size_t, 2>> edges(edge.nrow());
  for (size_t i = 0; i < edge.nrow(); i++) {
    std::array<size_t, 2> to_add = {static_cast<size_t>(edge(i, 0)),
                                    static_cast<size_t>(edge(i, 1))};
    edges[i] = to_add;
  }

  return calc_crown_age(edges, el);
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
std::string l_to_newick(const Rcpp::NumericMatrix& ltable_R,
                        bool drop_extinct) {
  auto ltable_cpp = convert_to_ltable(ltable_R);
  auto newick_string = ltable_to_newick(ltable_cpp, drop_extinct);
  return newick_string;
}

// old stuff
/*
 double calc_nltt_cpp_old(const Rcpp::NumericVector& brts_one,
 const Rcpp::NumericVector& brts_two) {

 try {
 std::vector<double> v1(brts_one.begin(), brts_one.end());
 std::vector<double> v2(brts_two.begin(), brts_two.end());

 auto nltt = calc_nltt(v1, v2);
 return nltt;
 } catch(std::exception &ex) {
 forward_exception_to_r(ex);
 } catch(...) {
 ::Rf_error("c++ exception (unknown reason)");
 }
 return NA_REAL;
 }
 */
