#ifndef COMPUTEMATRICES_H
#define COMPUTEMATRICES_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec computeYStar (arma::mat X, arma::vec curr_beta, arma::mat Z, arma::mat Dinv, arma::vec curr_u, arma::vec y);
arma::mat computeVmu (arma::vec mu, double r);
arma::mat computeW (arma::mat Dinv, arma::mat V);
arma::mat computeVStar (arma::mat Z, arma::mat G, arma::mat W);
arma::mat computePREML (arma::mat Vsinv, arma::mat X);
arma::mat initialiseG (Rcpp::List rlevels, arma::vec sigmas);

#endif
