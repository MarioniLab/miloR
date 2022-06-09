#ifndef COMPUTEMATRICES_H
#define COMPUTEMATRICES_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec computeYStar (arma::mat X, arma::vec curr_beta, arma::mat Z, arma::mat Dinv, arma::vec curr_u, arma::vec y);
arma::mat computeVmu (arma::vec mu, double r);
arma::mat computeW (double disp, arma::mat Dinv, arma::mat V);
arma::mat computeVStar (arma::mat Z, arma::mat G, arma::mat W);
arma::mat computePREML (arma::mat Vsinv, arma::mat X);
arma::mat initialiseG (Rcpp::List rlevels, arma::vec sigmas);
arma::mat initialiseG_G (Rcpp::List u_indices, arma::vec sigmas, arma::mat Kin);
arma::mat invGmat (Rcpp::List rlevels, arma::vec sigmas);
arma::mat invGmat_G (Rcpp::List u_indices, arma::vec sigmas, arma::mat Kin);
arma::mat subMatG (arma::vec u_index, double sigma, arma::mat broadcast);

#endif
