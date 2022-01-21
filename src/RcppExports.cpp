// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// branching_times
std::vector< double > branching_times(const Rcpp::List& phy);
RcppExport SEXP _treestats_branching_times(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(branching_times(phy));
    return rcpp_result_gen;
END_RCPP
}
// calc_beta_cpp
double calc_beta_cpp(const Rcpp::List& phy, double upper_lim, std::string algorithm, double abs_tol, double rel_tol);
RcppExport SEXP _treestats_calc_beta_cpp(SEXP phySEXP, SEXP upper_limSEXP, SEXP algorithmSEXP, SEXP abs_tolSEXP, SEXP rel_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy(phySEXP);
    Rcpp::traits::input_parameter< double >::type upper_lim(upper_limSEXP);
    Rcpp::traits::input_parameter< std::string >::type algorithm(algorithmSEXP);
    Rcpp::traits::input_parameter< double >::type abs_tol(abs_tolSEXP);
    Rcpp::traits::input_parameter< double >::type rel_tol(rel_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_beta_cpp(phy, upper_lim, algorithm, abs_tol, rel_tol));
    return rcpp_result_gen;
END_RCPP
}
// calc_colless_cpp
double calc_colless_cpp(const Rcpp::List phy, std::string normalization);
RcppExport SEXP _treestats_calc_colless_cpp(SEXP phySEXP, SEXP normalizationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type phy(phySEXP);
    Rcpp::traits::input_parameter< std::string >::type normalization(normalizationSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_colless_cpp(phy, normalization));
    return rcpp_result_gen;
END_RCPP
}
// calc_colless_cpp2
double calc_colless_cpp2(const std::vector< int >& edge, std::string normalization);
RcppExport SEXP _treestats_calc_colless_cpp2(SEXP edgeSEXP, SEXP normalizationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< int >& >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< std::string >::type normalization(normalizationSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_colless_cpp2(edge, normalization));
    return rcpp_result_gen;
END_RCPP
}
// calc_colless_cpp3
double calc_colless_cpp3(const Rcpp::NumericMatrix& edge, std::string normalization);
RcppExport SEXP _treestats_calc_colless_cpp3(SEXP edgeSEXP, SEXP normalizationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< std::string >::type normalization(normalizationSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_colless_cpp3(edge, normalization));
    return rcpp_result_gen;
END_RCPP
}
// calc_blum_cpp
double calc_blum_cpp(const Rcpp::List phy);
RcppExport SEXP _treestats_calc_blum_cpp(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(calc_blum_cpp(phy));
    return rcpp_result_gen;
END_RCPP
}
// calc_sackin_cpp
double calc_sackin_cpp(const Rcpp::List phy, const Rcpp::String& normalization);
RcppExport SEXP _treestats_calc_sackin_cpp(SEXP phySEXP, SEXP normalizationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type phy(phySEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type normalization(normalizationSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_sackin_cpp(phy, normalization));
    return rcpp_result_gen;
END_RCPP
}
// calc_nltt_cpp
double calc_nltt_cpp(const Rcpp::List& phy1, const Rcpp::List& phy2);
RcppExport SEXP _treestats_calc_nltt_cpp(SEXP phy1SEXP, SEXP phy2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy1(phy1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy2(phy2SEXP);
    rcpp_result_gen = Rcpp::wrap(calc_nltt_cpp(phy1, phy2));
    return rcpp_result_gen;
END_RCPP
}
// calc_gamma_cpp
double calc_gamma_cpp(const Rcpp::List& phy);
RcppExport SEXP _treestats_calc_gamma_cpp(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(calc_gamma_cpp(phy));
    return rcpp_result_gen;
END_RCPP
}
// calc_phylodiv_cpp
double calc_phylodiv_cpp(const Rcpp::List& phy, double t, double crown_age, double extinct_acc);
RcppExport SEXP _treestats_calc_phylodiv_cpp(SEXP phySEXP, SEXP tSEXP, SEXP crown_ageSEXP, SEXP extinct_accSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy(phySEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< double >::type extinct_acc(extinct_accSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_phylodiv_cpp(phy, t, crown_age, extinct_acc));
    return rcpp_result_gen;
END_RCPP
}
// calc_rho_cpp
double calc_rho_cpp(const Rcpp::List& phy, double crown_age);
RcppExport SEXP _treestats_calc_rho_cpp(SEXP phySEXP, SEXP crown_ageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy(phySEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_rho_cpp(phy, crown_age));
    return rcpp_result_gen;
END_RCPP
}
// phylo_to_l
Rcpp::NumericMatrix phylo_to_l(const Rcpp::List& phy);
RcppExport SEXP _treestats_phylo_to_l(SEXP phySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type phy(phySEXP);
    rcpp_result_gen = Rcpp::wrap(phylo_to_l(phy));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_treestats_branching_times", (DL_FUNC) &_treestats_branching_times, 1},
    {"_treestats_calc_beta_cpp", (DL_FUNC) &_treestats_calc_beta_cpp, 5},
    {"_treestats_calc_colless_cpp", (DL_FUNC) &_treestats_calc_colless_cpp, 2},
    {"_treestats_calc_colless_cpp2", (DL_FUNC) &_treestats_calc_colless_cpp2, 2},
    {"_treestats_calc_colless_cpp3", (DL_FUNC) &_treestats_calc_colless_cpp3, 2},
    {"_treestats_calc_blum_cpp", (DL_FUNC) &_treestats_calc_blum_cpp, 1},
    {"_treestats_calc_sackin_cpp", (DL_FUNC) &_treestats_calc_sackin_cpp, 2},
    {"_treestats_calc_nltt_cpp", (DL_FUNC) &_treestats_calc_nltt_cpp, 2},
    {"_treestats_calc_gamma_cpp", (DL_FUNC) &_treestats_calc_gamma_cpp, 1},
    {"_treestats_calc_phylodiv_cpp", (DL_FUNC) &_treestats_calc_phylodiv_cpp, 4},
    {"_treestats_calc_rho_cpp", (DL_FUNC) &_treestats_calc_rho_cpp, 2},
    {"_treestats_phylo_to_l", (DL_FUNC) &_treestats_phylo_to_l, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_treestats(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
