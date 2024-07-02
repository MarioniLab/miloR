#ifndef COMPUTEMATRICES_H
#define COMPUTEMATRICES_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec computeYStar (arma::mat X, arma::vec curr_beta, arma::mat Z, arma::mat Dinv,
                        arma::vec curr_u, arma::vec y, arma::vec offsets);
arma::mat computeVmu (arma::vec mu, double r, std::string vardist);
arma::mat computeVmuPoisson(arma::vec mu);
arma::mat computeVmuNB(arma::vec mu, double r);
arma::mat computeW (double disp, arma::mat Dinv, std::string vardist);
arma::mat computeWNB(double disp, arma::mat Dinv);
arma::mat computeWPoisson(arma::mat Dinv);
arma::mat computeVStar (arma::mat Z, arma::mat G, arma::mat W);
arma::mat computeBupdate(const arma::mat& Gdiff, const arma::mat& Z, const arma::mat& Wdiff);
arma::mat computePREML (const arma::mat& Vsinv, const arma::mat& X);
arma::mat initialiseG (Rcpp::List rlevels, arma::vec sigmas);
arma::mat initialiseG_G (Rcpp::List u_indices, arma::vec sigmas, arma::mat Kin);
arma::mat invGmat (Rcpp::List rlevels, arma::vec sigmas);
arma::mat invGmat_G (Rcpp::List u_indices, arma::vec sigmas, arma::mat Kin);
arma::mat subMatG (arma::vec u_index, double sigma, arma::mat broadcast);
// arma::mat makePCGFill(const Rcpp::List& u_indices, const arma::mat& Kinv);
arma::mat broadcastInverseMatrix(arma::mat matrix, const unsigned int& n);

#endif
