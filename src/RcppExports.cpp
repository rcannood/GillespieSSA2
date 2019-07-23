// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// make_ode_em
SEXP make_ode_em(double tau, double noise_strength);
RcppExport SEXP _gillespie_make_ode_em(SEXP tauSEXP, SEXP noise_strengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type noise_strength(noise_strengthSEXP);
    rcpp_result_gen = Rcpp::wrap(make_ode_em(tau, noise_strength));
    return rcpp_result_gen;
END_RCPP
}
// make_ssa_btl
SEXP make_ssa_btl(double mean_firings);
RcppExport SEXP _gillespie_make_ssa_btl(SEXP mean_firingsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mean_firings(mean_firingsSEXP);
    rcpp_result_gen = Rcpp::wrap(make_ssa_btl(mean_firings));
    return rcpp_result_gen;
END_RCPP
}
// make_ssa_etl
SEXP make_ssa_etl(double tau);
RcppExport SEXP _gillespie_make_ssa_etl(SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(make_ssa_etl(tau));
    return rcpp_result_gen;
END_RCPP
}
// make_ssa_exact
SEXP make_ssa_exact();
RcppExport SEXP _gillespie_make_ssa_exact() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(make_ssa_exact());
    return rcpp_result_gen;
END_RCPP
}
// test_ssa_step
List test_ssa_step(SEXP ssa_alg, NumericVector& state, NumericVector& propensity, IntegerVector& nu_i, IntegerVector& nu_p, IntegerVector& nu_x);
RcppExport SEXP _gillespie_test_ssa_step(SEXP ssa_algSEXP, SEXP stateSEXP, SEXP propensitySEXP, SEXP nu_iSEXP, SEXP nu_pSEXP, SEXP nu_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ssa_alg(ssa_algSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type state(stateSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type propensity(propensitySEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type nu_i(nu_iSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type nu_p(nu_pSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type nu_x(nu_xSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ssa_step(ssa_alg, state, propensity, nu_i, nu_p, nu_x));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_gillespie();

static const R_CallMethodDef CallEntries[] = {
    {"_gillespie_make_ode_em", (DL_FUNC) &_gillespie_make_ode_em, 2},
    {"_gillespie_make_ssa_btl", (DL_FUNC) &_gillespie_make_ssa_btl, 1},
    {"_gillespie_make_ssa_etl", (DL_FUNC) &_gillespie_make_ssa_etl, 1},
    {"_gillespie_make_ssa_exact", (DL_FUNC) &_gillespie_make_ssa_exact, 0},
    {"_gillespie_test_ssa_step", (DL_FUNC) &_gillespie_test_ssa_step, 6},
    {"_rcpp_module_boot_gillespie", (DL_FUNC) &_rcpp_module_boot_gillespie, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_gillespie(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
