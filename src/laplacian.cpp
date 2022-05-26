#include <Rcpp.h>
#include "util.h"
#include "dist_nodes.h"

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//' function to create laplacian matrix
//' @param phy phy
//' @return numericmatrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix prep_lapl_spec(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);

  std::vector< std::vector< double >> lapl_mat = dist_nodes(edge, el);
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
std::vector<double> laplacian_eigen_vals(const Rcpp::List& phy) {
  auto edge = phy_to_edge(phy);
  auto el   = phy_to_el(phy);

  std::vector< std::vector< double >> lapl_mat = dist_nodes(edge, el);

  for (size_t i = 0; i < lapl_mat.size(); ++i) {
    lapl_mat[i][i] = - std::accumulate(lapl_mat[i].begin(), lapl_mat[i].end(), 0.0);
  }

  Eigen::MatrixXd eigen_matrix(lapl_mat.size(),
                               lapl_mat[0].size());

  for (size_t i = 0; i < lapl_mat.size(); ++i) {
    for (size_t j = 0; j < lapl_mat[i].size(); ++j) {
      eigen_matrix(i, j) = lapl_mat[i][j];
    }
  }

  bool use_symmetric = true;

  if (use_symmetric) {

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

    solver.computeDirect(eigen_matrix, Eigen::ComputeEigenvectors);

    auto x = solver.eigenvalues();

    std::vector<double> ev;
    for (size_t i = 0; i < x.size(); ++i) {
      if (x[i] >= 1.0) ev.push_back(x[i]);
    }
    std::sort(ev.begin(), ev.end(), std::greater<double>());
    return ev;
  } else {

    Eigen::EigenSolver<Eigen::MatrixXd> solver;
    solver.compute(eigen_matrix, Eigen::ComputeEigenvectors);
    auto x = solver.eigenvalues();
    std::vector<double> ev;
    for (size_t i = 0; i < x.size(); ++i) {
      if (x[i].real() >= 1.0) ev.push_back(x[i].real());
    }
    std::sort(ev.begin(), ev.end(), std::greater<double>());
    return ev;
  }
}
