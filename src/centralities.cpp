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
  try {
    return diameter(edge, el, weight);
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









