// Network statistics based

#include <vector>
#include <array>
#include <Rcpp.h>

#include "centralities.h"
#include "util.h"

// [[Rcpp::export]]
double calc_wiener_cpp(const Rcpp::List& phy,
                       bool normalize,
                       bool weight) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return wiener(edge, el, normalize, weight);
}

// [[Rcpp::export]]
double calc_max_betweenness_cpp(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return max_betweenness(edge, el);
}

// [[Rcpp::export]]
double calc_max_betweenness_ltable_cpp(const Rcpp::NumericMatrix& l_from_R) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  return max_betweenness_ltable(l_in_cpp);
}

// [[Rcpp::export]]
double calc_max_closeness_cpp(const Rcpp::List& phy, bool weight) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return max_closeness(edge, el, weight);
}

// [[Rcpp::export]]
double calc_diameter_cpp(const Rcpp::List& phy, bool weight) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return diameter(edge, el, weight);
}

// [[Rcpp::export]]
double calc_diameter_ltable_cpp(const Rcpp::NumericMatrix& l_from_R,
                                bool weight) {
  auto l_in_cpp = convert_to_ltable(l_from_R);
  return diameter_ltable(l_in_cpp, weight);
}
