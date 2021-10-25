#include <Rcpp.h>

#include "beta.h"
#include "nltt.h"
#include "sackin.h"
#include "gamma.h"


// [[Rcpp::export]]
double calc_beta_cpp(Rcpp::NumericMatrix in_table,
                     double upper_lim) {
  if (in_table.ncol() != 4) {
    Rcpp::stop("ltable needs four columns");
  }
  ltable ltab;
  for (int i = 0; i < in_table.nrow(); ++i) {
    std::array< float, 4 > row;
    for (int j = 0; j < in_table.ncol(); ++j) {
      row[j] = in_table(i,j);
    }
    ltab.push_back(row);
  }

  //  force_output("collected into vector\n");

  try {
    double output = calc_beta(ltab, -2, upper_lim);
    return output;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// [[Rcpp::export]]
int calc_sackin_cpp(Rcpp::NumericMatrix in_table) {

  std::vector< std::vector< float >> ltab;
  for (int i = 0; i < in_table.nrow(); ++i) {
    std::vector< float > row;
    for (int j = 0; j < in_table.ncol(); ++j) {
      row.push_back(in_table(i,j));
    }
    ltab.push_back(row);
  }

  //  force_output("collected into vector\n");

  try {
    int output = calc_sackin(ltab);
    return output;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

// [[Rcpp::export]]
float calc_nltt_cpp(const Rcpp::NumericVector& brts_one,
                    const Rcpp::NumericVector& brts_two)
{
  std::vector<float> v1(brts_one.begin(), brts_one.end());
  std::vector<float> v2(brts_two.begin(), brts_two.end());

  auto nltt = calc_nltt(v1, v2);
  return nltt;
}

// [[Rcpp::export]]
float calc_gamma_cpp(const Rcpp::NumericVector& brts_in)
{
  std::vector<float> brts(brts_in.begin(), brts_in.end());
  gamma_stat gamma_stat_object(brts);

  return gamma_stat_object.calc_gamma_stat();
}
