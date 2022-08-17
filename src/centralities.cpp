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


// [[Rcpp::export]]
Rcpp::NumericMatrix get_adj_mat_cpp(const Rcpp::List& phy,
                                     bool weight) {
  auto edge = phy_to_edge(phy);
  Rcpp::NumericVector el   = phy["edge.length"];
  int num_nodes_tips = el.size() + 1;

  Rcpp::NumericMatrix out(num_nodes_tips, num_nodes_tips);

  if (weight) {
    size_t cnt = 0;
    for (const auto& i : edge) {
      int x = i[0] - 1;
      int y = i[1] - 1;
      out(x, y) = el(cnt);
      out(y, x) = el(cnt);
      cnt++;
    }

  } else {
    for (const auto& i : edge) {
      int x = i[0] - 1;
      int y = i[1] - 1;
      out(x, y) = 1.0;
      out(y, x) = 1.0;
    }
  }
  return out;
}









