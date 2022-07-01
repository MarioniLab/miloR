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
#endif
