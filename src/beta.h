#ifndef statistics_h
#define statistics_h

#include <random>
#include <cmath>
#include <array>
#include <unordered_map>
#include <nloptrAPI.h>



constexpr double pi_ = 3.14159265358979323846;

using ltable = std::vector< std::array<float, 4>>;

float max_map_size = 1e6;

float gammaln(float d)
{
  // return lgammaf(d);

  static std::unordered_map<float, float> cache;
  if (max_map_size == 0) return lgammaf(d);

  // look up function result
  float ret;
  auto it = cache.find(d);
  if (it != cache.end()) {
    ret = (*it).second;
  } else {
    ret = lgammaf(d);
    if (cache.size() < max_map_size) cache[d] = ret;
  }
  return ret;
}



class betastat {
public:
  betastat(const ltable& lt_in) : lt_(lt_in) {
    for (auto i : lt_) {
      brts_.push_back(i[0]);
    }
    std::sort(brts_.begin(), brts_.end());
    brts_.erase( std::unique(brts_.begin(), brts_.end()),
                 brts_.end());
    update_lr_matrix();
  };

  float calc_likelihood(float beta) const {
    std::vector< float > sn = get_sn(beta);
    std::vector< float > ll(lr_.size());

    for (size_t i = 0; i < n_.size(); ++i) {
      ll[i] = calc_log_prob(i, sn[ n_[i] ], beta);
    }
    double sumll = std::accumulate(ll.begin(), ll.end(), 0.f);
    return sumll;
  }

private:
  const ltable lt_;
  std::vector< std::array<size_t, 2>> lr_;
  std::vector< size_t > n_;
  std::vector< size_t > unique_n_;
  size_t max_n_;
  std::vector<float> brts_;

  int find_species_in_ltable(int sp) {
    for (int i = 0; i < static_cast<int>(lt_.size()); ++i) {
      if (lt_[i][2] == sp) {
        return i;
      }
    }
    return -1;
  }

  std::vector<float> find_daughters(int sp,
                                    float bt) {
    std::vector< float > output;
    for (auto i : lt_) {
      if (i[0] < bt && i[1] == sp) {
        output.push_back(i[2]);
      }
    }
    return(output);
  }

  size_t get_total_num_lin(int sp,
                           float bt) {

    int index = find_species_in_ltable(sp);
    size_t total_tips = 0;
    if (index >= 0) {
      if (lt_[index][3] == -1) {
        total_tips = 1;
      }
    }

    // now, we have to find daughters branching off
    std::vector< float > daughters = find_daughters(sp, bt);
    if (!daughters.empty())  {
      for (auto d : daughters) {
        total_tips += get_total_num_lin(static_cast<int>(d), bt);
      }
    }
    return(total_tips);
  }

  std::vector< size_t > get_indices(float bt) {
    std::vector< size_t > indices;
    for (size_t i = 0; i < lt_.size(); ++i) {
      if (lt_[i][0] == bt) {
        indices.push_back(i);
      }
    }
    return indices;
  }

  void update_lr_matrix() {
    max_n_ = 0;
    for (auto br : brts_) {
      std::array<size_t, 2> lr = {0, 0};
      std::vector< size_t > indices = get_indices(br);
      if (indices.size() == 2) {
        lr[0] = get_total_num_lin(lt_[indices[0]][2], br);
        lr[1] = get_total_num_lin(lt_[indices[1]][2], br);
      }
      if (indices.size() == 1) {

        lr[0] = get_total_num_lin(lt_[indices[0]][2], br);
        lr[1] = get_total_num_lin(lt_[indices[0]][1], br);
      }

      if (lr[0] > lr[1]) {
        std::swap(lr[0], lr[1]);
      }
      size_t total_num_lin = lr[0] + lr[1];
      if (total_num_lin > max_n_) max_n_ = total_num_lin;
      n_.push_back(total_num_lin);
      lr_.push_back(lr);
    }
    unique_n_ = n_;
    std::sort(unique_n_.begin(), unique_n_.end());
    unique_n_.erase( std::unique(unique_n_.begin(), unique_n_.end()),
                     unique_n_.end());

    return;
  }

  float calc_i_n_b(size_t i, size_t n, float b) const {
    float nom = std::tgamma(1.f*(i + 1 + b)) * std::tgamma(1.f*(n - i + 1 + b));
    float denom = std::tgamma(1.f*(i + 1)) * std::tgamma(1.f*(n - i + 1));
    return(nom / denom);
  }

  float calc_i_n_b_l(size_t i, size_t n, float b) const {
    return gammaln(i + 1 + b) + gammaln(n - i + 1 + b) -
      gammaln(i + 1) - gammaln(n - i + 1);
  }



  std::vector<float> get_sn(float b) const {
    std::vector<float> sn(max_n_ + 1, 0.f);

    for (const auto& n : unique_n_) {
      for (size_t i = 1; i <= n - 1; ++i) {
        sn[n] += expf(calc_i_n_b_l(i, n, b));
      }
    }

    return sn;
  }

  float calc_log_prob(size_t index, float sn, float beta) const {
    float l = lr_[index][0];
    float r = lr_[index][1];
    return gammaln(beta + l + 1) + gammaln(beta + r + 1) -
      gammaln(l + 1) - gammaln(r + 1) - log(sn);
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


double calc_beta(const ltable& ltab,
                 double lower_lim,
                 double upper_lim) {
  betastat beta_calc(ltab);
  // now we do optimization?

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

