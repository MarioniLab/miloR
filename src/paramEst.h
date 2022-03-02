#ifndef PARAMEST_H
#define PARAMEST_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec sigmaScoreREML (Rcpp::List pvstar_i, arma::mat Vsinv, arma::vec ystar, arma::mat P);
arma::mat sigmaInfoREML (Rcpp::List pvstari);
arma::vec sigmaScore (arma::vec ystar, arma::vec beta, arma::mat X, Rcpp::List V_partial, arma::mat V_star_inv);
arma::mat sigmaInformation (arma::mat V_star_inv, Rcpp::List V_partial);
arma::vec FisherScore (arma::mat hess, arma::vec score_vec, arma::vec theta_hat);
arma::vec solve_equations (arma::mat X, arma::mat Winv, arma::mat Z, arma::mat Ginv, arma::vec beta, arma::vec u, arma::vec ystar);

#endif
