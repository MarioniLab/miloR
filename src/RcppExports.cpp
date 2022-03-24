// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// computeVStar
arma::mat computeVStar(arma::mat Z, arma::mat G, arma::mat W);
RcppExport SEXP _miloR_computeVStar(SEXP ZSEXP, SEXP GSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(computeVStar(Z, G, W));
    return rcpp_result_gen;
END_RCPP
}
// fitPLGlmm
List fitPLGlmm(const arma::mat& Z, const arma::mat& X, arma::vec muvec, arma::vec curr_beta, arma::vec curr_theta, arma::vec curr_u, arma::vec curr_sigma, arma::mat curr_G, const arma::vec& y, List u_indices, double theta_conv, const List& rlevels, double curr_disp, const bool& REML, const int& maxit);
RcppExport SEXP _miloR_fitPLGlmm(SEXP ZSEXP, SEXP XSEXP, SEXP muvecSEXP, SEXP curr_betaSEXP, SEXP curr_thetaSEXP, SEXP curr_uSEXP, SEXP curr_sigmaSEXP, SEXP curr_GSEXP, SEXP ySEXP, SEXP u_indicesSEXP, SEXP theta_convSEXP, SEXP rlevelsSEXP, SEXP curr_dispSEXP, SEXP REMLSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muvec(muvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type curr_beta(curr_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type curr_theta(curr_thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type curr_u(curr_uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type curr_sigma(curr_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type curr_G(curr_GSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type u_indices(u_indicesSEXP);
    Rcpp::traits::input_parameter< double >::type theta_conv(theta_convSEXP);
    Rcpp::traits::input_parameter< const List& >::type rlevels(rlevelsSEXP);
    Rcpp::traits::input_parameter< double >::type curr_disp(curr_dispSEXP);
    Rcpp::traits::input_parameter< const bool& >::type REML(REMLSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(fitPLGlmm(Z, X, muvec, curr_beta, curr_theta, curr_u, curr_sigma, curr_G, y, u_indices, theta_conv, rlevels, curr_disp, REML, maxit));
    return rcpp_result_gen;
END_RCPP
}
// invertPseudoVar
arma::mat invertPseudoVar(arma::mat A, arma::mat B, arma::mat Z);
RcppExport SEXP _miloR_invertPseudoVar(SEXP ASEXP, SEXP BSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(invertPseudoVar(A, B, Z));
    return rcpp_result_gen;
END_RCPP
}
// multiP
List multiP(List partials, arma::mat psvar_in);
RcppExport SEXP _miloR_multiP(SEXP partialsSEXP, SEXP psvar_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type partials(partialsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type psvar_in(psvar_inSEXP);
    rcpp_result_gen = Rcpp::wrap(multiP(partials, psvar_in));
    return rcpp_result_gen;
END_RCPP
}
// pseudovarPartial
List pseudovarPartial(arma::mat x, List rlevels, StringVector cnames);
RcppExport SEXP _miloR_pseudovarPartial(SEXP xSEXP, SEXP rlevelsSEXP, SEXP cnamesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type rlevels(rlevelsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type cnames(cnamesSEXP);
    rcpp_result_gen = Rcpp::wrap(pseudovarPartial(x, rlevels, cnames));
    return rcpp_result_gen;
END_RCPP
}
// pseudovarPartial_C
List pseudovarPartial_C(arma::mat Z, List u_indices);
RcppExport SEXP _miloR_pseudovarPartial_C(SEXP ZSEXP, SEXP u_indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< List >::type u_indices(u_indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(pseudovarPartial_C(Z, u_indices));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_miloR_computeVStar", (DL_FUNC) &_miloR_computeVStar, 3},
    {"_miloR_fitPLGlmm", (DL_FUNC) &_miloR_fitPLGlmm, 15},
    {"_miloR_invertPseudoVar", (DL_FUNC) &_miloR_invertPseudoVar, 3},
    {"_miloR_multiP", (DL_FUNC) &_miloR_multiP, 2},
    {"_miloR_pseudovarPartial", (DL_FUNC) &_miloR_pseudovarPartial, 3},
    {"_miloR_pseudovarPartial_C", (DL_FUNC) &_miloR_pseudovarPartial_C, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_miloR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}