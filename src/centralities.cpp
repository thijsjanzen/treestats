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
double calc_max_closeness_cpp(const Rcpp::List& phy, bool weight) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  try {
    return max_closeness(edge, el, weight);
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
double calc_diameter_cpp(const Rcpp::List& phy, bool weight) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);
  return diameter(edge, el, weight);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix get_adj_mat_cpp(const std::vector<long>& parent_list,
                                    const std::vector<double>& bl,
                                    bool weight) {

  int num_nodes_tips = *std::max_element(parent_list.begin(), parent_list.end());

  Rcpp::NumericMatrix out(num_nodes_tips, num_nodes_tips);

  for (size_t i = 0; i < parent_list.size(); i += 2 ) {

    int x = static_cast<int>(parent_list[i]) - 1;
    int y = static_cast<int>(parent_list[i + 1]) - 1;

    double val = weight ? bl[i / 2] : 1.0;
    out(x, y) = val;
    out(y, x) = val;
  }
  return out;
}



