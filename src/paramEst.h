#ifndef PARAMEST_H
#define PARAMEST_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

arma::vec sigmaScoreREML_arma (const Rcpp::List& pvstar_i, const arma::vec& ystar,
                               const arma::mat& P, const arma::vec& curr_beta,
                               const arma::mat& X, const arma::mat& Vstarinv,
                               const Rcpp::List& remldiffV);
arma::mat sigmaInfoREML_arma (const Rcpp::List& pvstari, const arma::mat& P);
arma::vec sigmaScore (arma::vec ystar, arma::vec beta, arma::mat X, Rcpp::List V_partial, arma::mat V_star_inv);
arma::mat sigmaInformation (arma::mat V_star_inv, Rcpp::List V_partial);
arma::vec fisherScore (const arma::mat& hess, const arma::vec& score_vec, const arma::vec& theta_hat);
arma::vec solveEquations (const int& c, const int& m, const arma::mat& ZtWinv, const arma::mat& XtWinv,
                          const arma::mat& coeffmat, const arma::vec& beta, const arma::vec& u, const arma::vec& ystar);
// arma::vec solveEquationsPCG (const int& c, const int& m, const arma::mat& Winv, const arma::mat& Zt, const arma::mat& Xt,
//                              const arma::mat& coeffmat, const arma::vec& curr_theta, const arma::vec& ystar, const double& conv_tol);
arma::mat coeffMatrix(const arma::mat& X, const arma::mat& XtWinv, const arma::mat& ZtWinv,
                      const arma::mat& Z, const arma::mat& Ginv);
arma::mat computeZstar(const arma::mat& Z, const arma::vec& curr_sigma, const Rcpp::List& u_indices);
// arma::vec conjugateGradient(const arma::mat& A, const arma::vec& x, const arma::vec& b, double conv_tol);
arma::vec estHasemanElston(const arma::mat& Z, const arma::mat& PREML,
                           const Rcpp::List& u_indices, const arma::vec& ystar,
                           const arma::mat& PZ);
arma::vec estHasemanElstonML(const arma::mat& Z, const Rcpp::List& u_indices,
                             const arma::vec& ystar);
arma::vec estHasemanElstonGenetic(const arma::mat& Z, const arma::mat& PREML, const arma::mat& PZ,
                                  const Rcpp::List& u_indices, const arma::vec& ystar,
                                  const arma::mat& Kin);
arma::vec estHasemanElstonConstrained(const arma::mat& Z, const arma::mat& PREML, const Rcpp::List& u_indices,
                                      const arma::vec& ystar, arma::vec he_update, const int& Iters,
                                      const arma::mat& PZ);
arma::vec estHasemanElstonConstrainedML(const arma::mat& Z, const Rcpp::List& u_indices,
                                        const arma::vec& ystar, arma::vec he_update, const int& Iters);
arma::vec estHasemanElstonConstrainedGenetic(const arma::mat& Z, const arma::mat& PREML, const arma::mat& PZ,
                                             const Rcpp::List& u_indices,
                                             const arma::vec& ystar, const arma::mat& Kin, arma::vec he_update,
                                             const int& Iters);
arma::vec nnlsSolve(const arma::mat& vecZ, const arma::vec& Y, arma::vec nnls_update, const int& Iters);
arma::vec fastNnlsSolve(const arma::mat& vecZ, const arma::vec& Y);
arma::mat vectoriseZ(const arma::mat& Z, const Rcpp::List& u_indices, const arma::mat& P,
                     const arma::mat& PZ);
arma::mat vectoriseZML(const arma::mat& Z, const Rcpp::List& u_indices);
arma::mat vectoriseZGenetic(const arma::mat& Z, const Rcpp::List& u_indices,
                            const arma::mat& P, const arma::mat& PZ, const arma::mat& Kin);
double phiLineSearch(double disp, double lower, double upper, const int& c,
                     const arma::vec& mu, const arma::mat& Ginv, double pi,
                     const arma::vec& curr_u, const arma::vec& sigma,
                     const arma::vec& y);
double phiGoldenSearch(double disp, double lower, double upper, const int& c,
                       const arma::vec& mu, const arma::mat& Ginv, double pi,
                       const arma::vec& curr_u, const arma::vec& sigma,
                       const arma::vec& y);
double phiMME(const arma::vec& y, const arma::vec& curr_sigma);
double nbLogLik(const arma::vec& mu, double phi, const arma::vec& y);
double normLogLik(const int& c, const arma::mat& Ginv, const arma::mat& G,
                  const arma::vec& curr_u, double pi);
#endif
