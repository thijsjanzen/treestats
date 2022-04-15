// DEPTH BASED

#include <vector>
#include <array>
#include <Rcpp.h>

#include "util.h"
#include "depth.h"
#include "ltable.h"

// [[Rcpp::export]]
int calc_max_del_width_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker tree(parent_list);
  return tree.calc_max_del_width();
}

// [[Rcpp::export]]
double calc_max_del_width_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.max_del_width();
}



// [[Rcpp::export]]
int calc_max_width_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker tree(parent_list);
  return tree.calc_max_width();
}

// [[Rcpp::export]]
double calc_max_width_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.calc_max_width();
}

// [[Rcpp::export]]
int calc_max_depth_cpp(const std::vector<long>& parent_list) {
  depth::phylo_tree local_tree(parent_list);
  return local_tree.max_depth();
}

// [[Rcpp::export]]
double calc_max_depth_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat s(l_in_cpp);
  return s.calc_max_depth();
}


// [[Rcpp::export]]
double calc_var_leaf_depth_cpp(const std::vector<long>& parent_list) {
  width::depth_tracker local_tree(parent_list);
  return local_tree.var_leaf_depth();
}

// [[Rcpp::export]]
double calc_var_leaf_depth_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  ltab::stat local_tree(l_in_cpp);
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
  ltab::stat s(l_in_cpp);
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
  ltab::stat s(l_in_cpp);
  return s.calc_b2();
}


