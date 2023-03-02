#ifndef PARAMEST_H
#define PARAMEST_H

#include<RcppArmadillo.h>
#include<RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

arma::vec sigmaScoreREML_arma (Rcpp::List pvstar_i, const arma::vec& ystar, const arma::mat& P);
arma::mat sigmaInfoREML_arma (const Rcpp::List& pvstari, const arma::mat& P);
arma::vec sigmaScore (arma::vec ystar, arma::vec beta, arma::mat X, Rcpp::List V_partial, arma::mat V_star_inv);
arma::mat sigmaInformation (arma::mat V_star_inv, Rcpp::List V_partial);
arma::vec fisherScore (arma::mat hess, arma::vec score_vec, arma::vec theta_hat);
arma::vec solveEquations (const int& c, const int& m, const arma::mat& Winv, const arma::mat& Zt, const arma::mat& Xt,
                          arma::mat coeffmat, arma::vec beta, arma::vec u, const arma::vec& ystar);
arma::vec solveEquationsPCG (const int& c, const int& m, const arma::mat& Winv, const arma::mat& Zt, const arma::mat& Xt,
                             arma::mat coeffmat, arma::vec curr_theta, const arma::vec& ystar, const double& conv_tol);
arma::mat coeffMatrix(const arma::mat& X, const arma::mat& Winv, const arma::mat& Z, const arma::mat& Ginv);
arma::mat computeZstar(const arma::mat& Z, arma::vec curr_sigma, Rcpp::List u_indices);
arma::vec conjugateGradient(arma::mat A, arma::vec x, arma::vec b, double conv_tol);
arma::vec estHasemanElston(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices, arma::vec ystar);
arma::vec estHasemanElstonGenetic(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices, arma::vec ystar, arma::mat Kin);
arma::vec estHasemanElstonConstrained(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices,
                                      arma::vec ystar, arma::vec he_update, const int& Iters);
arma::vec estHasemanElstonConstrainedGenetic(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices,
                                             arma::vec ystar, arma::mat Kin, arma::vec he_update, const int& Iters);
arma::vec nnlsSolve(const arma::mat& vecZ, arma::vec Y, arma::vec init_nnls, const int& Iters);
arma::vec fastNnlsSolve(const arma::mat& vecZ, arma::vec Y);
arma::mat vectoriseZ(arma::mat Z, Rcpp::List u_indices, arma::mat P);
arma::mat vectoriseZGenetic(arma::mat Z, Rcpp::List u_indices, arma::mat P, arma::mat Kin);
double phiLineSearch(double disp, double lower, double upper, const int& c,
                     arma::vec mu, arma::mat Ginv, double pi,
                     arma::vec curr_u, arma::vec sigma, arma::vec y);
double phiGoldenSearch(double disp, double lower, double upper, const int& c,
                       arma::vec mu, arma::mat Ginv, double pi,
                       arma::vec curr_u, arma::vec sigma, arma::vec y);
double nbLogLik(arma::vec mu, double phi, arma::vec y);
double normLogLik(const int& c, arma::mat Ginv, arma::mat G, arma::vec curr_u, double pi);
#endif
