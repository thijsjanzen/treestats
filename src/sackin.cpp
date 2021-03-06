// SACKIN BASED

#include <vector>
#include <array>
#include <Rcpp.h>

#include "util.h"
#include "sackin.h"
#include "ltable.h"
#include "cherries.h"

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
  return calc_sackin(local_ltab, normalization);
}


// [[Rcpp::export]]
double calc_tot_coph_cpp(const std::vector<long>& tree_edge) {
  phylo_tree sackin_tree(tree_edge);
  return sackin_tree.calc_tot_coph();
}

// [[Rcpp::export]]
double calc_tot_coph_ltable_cpp(const Rcpp::NumericMatrix& ltab) {
  auto local_ltab = convert_to_ltable(ltab);
  ltab::stat s(local_ltab);
  return s.calc_tot_coph();
}

// [[Rcpp::export]]
double calc_blum_cpp(const std::vector<long>& tree_edge) {
  phylo_tree sackin_tree(tree_edge);
  return sackin_tree.calc_blum();
}
// [[Rcpp::export]]
double calc_blum_ltable_cpp(const Rcpp::NumericMatrix& ltab_in) {
  auto local_ltab = convert_to_ltable(ltab_in);
  ltab::stat s(local_ltab);
  return s.calc_blum();
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
size_t pitchforks_cpp(const std::vector<long>& tree_edge) {
  // ltable version uses colless
  phylo_tree sackin_tree(tree_edge);
  return sackin_tree.count_pitchforks();
}

