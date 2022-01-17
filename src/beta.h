#ifndef beta_h
#define beta_h

#include <random>
#include <cmath>
#include <array>
#include <unordered_map>
#include <nloptrAPI.h>
#include <algorithm>


class betastat {
public:
  betastat(const std::vector< std::array< size_t, 2 >> e) : edge(e) {
    tiplist = std::vector<int>(edge.size() + 2, -1);
    update_lr_matrix();
  }

  double calc_likelihood(double beta) const {
    std::vector< double > sn = get_sn(beta);
    std::vector< double > ll(lr_.size());
    assert(ll.size() == n_.size());
    for (size_t i = 0; i < n_.size(); ++i) {
      auto index = n_[i];
      assert(index < sn.size() );
      ll[i] = calc_log_prob(i, sn[ index ], beta);
    }
    double sumll = std::accumulate(ll.begin(), ll.end(), 0.f);
    return sumll;
  }

private:
  std::vector< std::array<size_t, 2>> lr_;
  std::vector< std::array<size_t, 2>> edge;
  size_t max_n_;
  std::vector< size_t > n_;

  std::vector< int > tiplist;

  size_t get_num_tips(size_t label, size_t root_label) {
    if (label >= tiplist.size()) {
      throw std::out_of_range("label > tiplist.size()");
    }

    if (label < root_label) {
      return 1;
    }

    if (tiplist[label] > 0) {
      return(tiplist[label]);
    }

    std::vector< size_t > matches;
    auto match1 = std::lower_bound(edge.begin(), edge.end(), label, [&](const auto& a, size_t val){
      return a[0] < val;
    });

    while (match1 != edge.end()) {
      if ((*match1)[0] == label) {
        matches.push_back((*match1)[1]);
        match1++;
      } else {
        break;
      }
    }

    if (matches.empty()) {
      // this can't really happen.
      tiplist[label] = 1;
      return 1;
    }

    size_t s = 0;
    for (size_t j = 0; j < matches.size(); ++j) {
      s += get_num_tips(matches[j], root_label);
    }
    tiplist[label] = s;
    return s;
  }


  void update_lr_matrix() {

    size_t root_label = edge[0][0];

    std::sort(edge.begin(), edge.end(), [&](const auto& a, const auto& b) {
      return a[0] < b[0];
    });


    for (size_t i = 0; i < edge.size(); ++i) {
      if (i + 1 < edge.size()) {
        auto j = i + 1;
        std::array< size_t, 2> lr;
        lr[0] = get_num_tips(edge[i][1], root_label);
        lr[1] = get_num_tips(edge[j][1], root_label);

        if (lr[0] > lr[1]) {
          //std::swap(lr[0], lr[1]);
          size_t temp = lr[1];
          lr[1] = lr[0];
          lr[0] = temp;
        }
        size_t total_num_lin = lr[0] + lr[1];
        n_.push_back(total_num_lin);
        lr_.push_back(lr);
        i++;
      }
    }

    max_n_ = *std::max_element(n_.cbegin(), n_.cend());
  }


  double calc_i_n_b(size_t i, size_t n, double b) const {
    double nom = std::tgamma(1.f*(i + 1 + b)) * std::tgamma(1.f*(n - i + 1 + b));
    double denom = std::tgamma(1.f*(i + 1)) * std::tgamma(1.f*(n - i + 1));
    return(nom / denom);
  }

  double calc_i_n_b_l(size_t i, size_t n, double b) const {
    return lgamma(i + 1 + b) + lgamma(n - i + 1 + b) -
      lgamma(i + 1) - lgamma(n - i + 1);
  }

  std::vector<double> get_sn(double b) const {
    std::vector<double> sn(max_n_ + 1, 0.0);
    std::vector<double> xn(max_n_ + 1, 0.0);

    xn[2] = 1.0;
    xn[3] = 0.5;

    sn[2] = expf(calc_i_n_b_l(1, 2, b));
    sn[3] = expf(calc_i_n_b_l(1, 3, b)) + expf(calc_i_n_b_l(2, 3, b));

    for (size_t n = 3; n < max_n_; ++n) {

      auto term1 = n + 2 + 2 * b;
      auto term2 = 2 * (n + b) * xn[n];

      xn[n + 1] = ((n + b) * (n + 1) * xn[n]) / (n * term1 + term2);
      sn[n + 1] =  (1.0 / (n + 1) ) * (term1 + term2 / n) * sn[n];
    }

    return sn;
  }

  double calc_log_prob(size_t index, double sn, double beta) const {
    double l = lr_[index][0];
    double r = lr_[index][1];
    return lgamma(beta + l + 1) + lgamma(beta + r + 1) -
           lgamma(l + 1) - lgamma(r + 1) - log(sn);
  }
};

struct nlopt_f_data {

  nlopt_f_data(const betastat& b_in) : b(b_in) {
  }

  const betastat b;
};

double objective(unsigned int n, const double* x, double*, void* func_data) {
  auto psd = reinterpret_cast<nlopt_f_data*>(func_data);
  return(-psd->b.calc_likelihood(x[0]));
}


double calc_beta(const std::vector< std::array< size_t, 2 >>& edge,
                 double lower_lim,
                 double upper_lim) {
  betastat beta_calc(edge);
  // now we do optimization

  nlopt_f_data optim_data(beta_calc);

  nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX, static_cast<unsigned int>(1));
  double llim[1] = {static_cast<double>(lower_lim)};
  double ulim[1] = {static_cast<double>(upper_lim)};

  nlopt_set_lower_bounds(opt, llim);
  nlopt_set_upper_bounds(opt, ulim);

  nlopt_set_min_objective(opt, objective, &optim_data);

  nlopt_set_xtol_rel(opt, 1e-6);
  std::vector<double> x = {0};
  double minf;

  auto nloptresult = nlopt_optimize(opt, &(x[0]), &minf);

  if (nloptresult < 0) {
    Rcpp::Rcout << "failure to optimize!\n";
  }

  nlopt_destroy(opt);

  double beta = x[0];
 // double ll = minf;
  return beta;
}




#endif /* statistics_h */

