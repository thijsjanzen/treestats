#include <Rcpp.h>

#include <cstring>


#include "branching_times.h"
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
#include "cherries.h"
#include "ILnumber.h"
#include "newick_to_edge.h"
#include "avgladder.h"
#include "depth.h"
#include "mntd.h"

using ltable = std::vector< std::array<double, 4>>;
using edge_table = std::vector< std::array< size_t, 2 >>;

auto convert_to_ltable(const Rcpp::NumericMatrix& mat_in) {
  ltable out(mat_in.nrow());

  for (size_t i = 0; i < mat_in.nrow(); ++i) {
    std::array<double, 4> row_entry = {mat_in(i, 0), mat_in(i, 1),
                                       mat_in(i, 2), mat_in(i, 3) };
    out[i] = row_entry;
  }
  return out;
}

// [[Rcpp::export]]
std::vector<double> branching_times_ltable_cpp(const Rcpp::NumericMatrix& mat_in) {
  std::vector<double> out(mat_in.nrow() - 1);
  for (size_t i = 1; i < mat_in.nrow(); ++i) {
    out[i - 1] = mat_in(i, 0);
  }
  return out;
}

// [[Rcpp::export]]
std::vector< double > branching_times_cpp(const Rcpp::List& phy) {

  std::vector< double > edge_length = phy["edge.length"];
  Rcpp::NumericMatrix edge = phy["edge"];

  size_t Nnode = phy["Nnode"];

  edge_table edge_cpp(edge.nrow());
  for (size_t i = 0; i < edge.nrow(); ++i) {
    std::array< size_t, 2 > row_entry = {static_cast<size_t>(edge(i, 0)),
                                         static_cast<size_t>(edge(i, 1))};
    edge_cpp[i] = row_entry;
  }

  sort_edge_and_edgelength(edge_cpp, edge_length);

  return branching_times(edge_cpp, edge_length, Nnode);
}

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

// [[Rcpp::export]]
double calc_colless_cpp(const std::vector<long>& parent_list,
                        std::string normalization) {
  colless_tree::phylo_tree colless_tree(parent_list);
  double output = static_cast<double>(colless_tree.calc_colless());

  if (normalization == "yule") {
    size_t n  = parent_list.size() / 4 + 1;
    output = colless_tree.correct_yule(output, n);
  }
  if (normalization == "pda") {
    size_t n  = parent_list.size() / 4 + 1;
    output = colless_tree.correct_pda(output, n);
  }
  return output;
}


// [[Rcpp::export]]
double calc_eWcolless_cpp(const std::vector<long>& parent_list) {
  colless_tree::phylo_tree colless_tree(parent_list);
  return colless_tree.calc_eWcolless();
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
double calc_Ibased_cpp(const std::vector<long>& parent_list) {
  colless_tree::phylo_tree colless_tree(parent_list);
  auto I_vec = colless_tree.collect_I();
  auto sum = std::accumulate(I_vec.begin(), I_vec.end(), 0.0);
  return sum * 1.0 / I_vec.size();
}

// [[Rcpp::export]]
double calc_Ibased_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  auto I_vec = s.collect_I();
  auto sum = std::accumulate(I_vec.begin(), I_vec.end(), 0.0);
  return sum * 1.0 / I_vec.size();
}



// [[Rcpp::export]]
int calc_rogers_cpp(const std::vector<long>& parent_list) {
  colless_tree::phylo_tree colless_tree(parent_list);
  return colless_tree.calc_rogers();
}

// [[Rcpp::export]]
double calc_rogers_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_rogers();
}

// [[Rcpp::export]]
double calc_var_leaf_depth_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker local_tree(parent_list);
  return local_tree.var_leaf_depth();
}

// [[Rcpp::export]]
double calc_var_leaf_depth_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  depth::depth_stat_ltab local_tree(l_in_cpp);
  return local_tree.calc_var_leaf_depth();
}

// [[Rcpp::export]]
int calc_sym_nodes_cpp(const std::vector<long>& parent_list) {

  width::depth_tracker tree(parent_list);
  return tree.calc_sym_nodes();
}

// [[Rcpp::export]]
double calc_b1_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker tree(parent_list);
  return tree.calc_b1();
}

// [[Rcpp::export]]
double calc_b1_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  depth::depth_stat_ltab s(l_in_cpp);
  return s.calc_b1();
}

// [[Rcpp::export]]
double calc_b2_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker tree(parent_list);
  return tree.calc_b2();
}

// [[Rcpp::export]]
double calc_b2_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  depth::depth_stat_ltab s(l_in_cpp);
  return s.calc_b2();
}



// [[Rcpp::export]]
int calc_max_depth_cpp(const std::vector<long>& parent_list) {
  depth::phylo_tree local_tree(parent_list);
  return local_tree.max_depth();
}

// [[Rcpp::export]]
double calc_max_depth_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  depth::depth_stat_ltab s(l_in_cpp);
  return s.calc_max_depth();
}

// [[Rcpp::export]]
int calc_max_width_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker tree(parent_list);
  return tree.calc_max_width();
}

// [[Rcpp::export]]
double calc_max_width_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  depth::depth_stat_ltab s(l_in_cpp);
  return s.calc_max_width();
}

// [[Rcpp::export]]
int calc_max_del_width_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker tree(parent_list);
  return tree.calc_max_del_width();
}

// [[Rcpp::export]]
double calc_max_del_width_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  depth::depth_stat_ltab s(l_in_cpp);
  return s.max_del_width();
}


// [[Rcpp::export]]
double calc_eWcolless_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  colless_stat_ltable s(l_in_cpp);
  return s.calc_ew_colless();
}

// [[Rcpp::export]]
double calc_blum_ltable_cpp(const Rcpp::NumericMatrix& ltab_in) {
  auto local_ltab = convert_to_ltable(ltab_in);
  sackin_stat_ltab s(local_ltab);
  return s.calc_blum();
}


// [[Rcpp::export]]
double calc_blum_cpp(const std::vector<long>& tree_edge) {
  phylo_tree sackin_tree(tree_edge);
  return sackin_tree.calc_blum();
}

// [[Rcpp::export]]
double calc_tot_coph_cpp(const std::vector<long>& tree_edge) {
  phylo_tree sackin_tree(tree_edge);
  return sackin_tree.calc_tot_coph();
}

// [[Rcpp::export]]
double calc_tot_coph_ltable_cpp(const Rcpp::NumericMatrix& ltab) {
  auto local_ltab = convert_to_ltable(ltab);
  sackin_stat_ltab s(local_ltab);
  return s.calc_tot_coph();
}


// [[Rcpp::export]]
double calc_sackin_cpp(const std::vector<long>& tree_edge,
                       const Rcpp::String& normalization) {

  phylo_tree sackin_tree(tree_edge);
  double output = static_cast<double>(sackin_tree.calc_sackin());

  if (normalization == "yule") {
    size_t n  = tree_edge.size() / 4 + 1;
    output = sackin_tree.correct_yule(n, output);
  }
  if (normalization == "pda") {
    size_t n  = tree_edge.size() / 4 + 1;
    output = sackin_tree.correct_pda(n, output);
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
  auto brts_one = branching_times_ltable_cpp(ltab1);
  auto brts_two = branching_times_ltable_cpp(ltab2);
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

  std::vector<double> brts_one = branching_times_cpp(phy1);
  std::vector<double> brts_two = branching_times_cpp(phy2);
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
  std::vector<double> brts = branching_times_cpp(phy);
  return calc_gamma(brts);
}

// [[Rcpp::export]]
double calc_gamma_ltable_cpp(const Rcpp::NumericMatrix& ltab_in) {
  std::vector<double> brts = branching_times_ltable_cpp(ltab_in);
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
    edge_table edges(edge.nrow());
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
  edge_table edges(edge.nrow());
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

  auto brts = branching_times_cpp(phy);
  return calc_rho(brts);
}

// [[Rcpp::export]]
double calc_rho_ltable_cpp(const Rcpp::NumericMatrix& ltab) {
  auto brts = branching_times_ltable_cpp(ltab);
  return calc_rho(brts);
}


// [[Rcpp::export]]
double calc_crown_age_cpp(const Rcpp::List& phy) {
  Rcpp::NumericMatrix edge = phy["edge"];
  Rcpp::NumericVector edge_length = phy["edge.length"];

  std::vector<double> el(edge_length.begin(), edge_length.end());
  edge_table edges(edge.nrow());
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
double calc_mntd_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto ltable = convert_to_ltable(ltable_R);
  return calc_mntd_ltable(ltable);
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
std::string l_to_newick(const Rcpp::NumericMatrix& ltable_R,
                        bool drop_extinct) {
  auto ltable_cpp = convert_to_ltable(ltable_R);
  auto newick_string = ltable_to_newick(ltable_cpp, drop_extinct);
  return newick_string;
}

// [[Rcpp::export]]
size_t pitchforks_cpp(const std::vector<long>& tree_edge) {
  phylo_tree sackin_tree(tree_edge);
  return sackin_tree.count_pitchforks();
}

// [[Rcpp::export]]
size_t pitchforks_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_pitchforks();
}

// [[Rcpp::export]]
size_t cherries_cpp(const std::vector<long>& tree_edge) {
  phylo_tree sackin_tree(tree_edge);
  return sackin_tree.count_cherries();
}

// [[Rcpp::export]]
size_t cherries_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  return calc_cherries_ltable(local_ltab);
}

// [[Rcpp::export]]
size_t ILnumber_cpp(const std::vector<long>& tree_edge) {
  return calc_IL(tree_edge);
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
double avgLadder_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto newick_str = l_to_newick(ltable_R, false);
  auto edge_table = newick_to_edge(newick_str);
  return avgLadder_cpp(edge_table);
}

// [[Rcpp::export]]
size_t ILnumber_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_IL();
}

// [[Rcpp::export]]
double stairs_cpp(const std::vector<long>& tree_edge) {
  colless_tree::phylo_tree tree(tree_edge);
  return tree.calc_stairs();
}

// [[Rcpp::export]]
double stairs2_cpp(const std::vector<long>& tree_edge) {
  colless_tree::phylo_tree tree(tree_edge);
  return tree.calc_stairs2();
}


// [[Rcpp::export]]
double stairs_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_stairs();
}

// [[Rcpp::export]]
double stairs2_ltable_cpp(const Rcpp::NumericMatrix& ltable_R) {
  auto local_ltab = convert_to_ltable(ltable_R);
  colless_stat_ltable c(local_ltab);
  return c.count_stairs2();
}
