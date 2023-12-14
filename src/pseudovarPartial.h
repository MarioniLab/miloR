#ifndef PSEUDOVARPARTIAL_H
#define PSEUDOVARPARTIAL_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List pseudovarPartial(arma::mat x, Rcpp::List rlevels, Rcpp::StringVector cnames);
Rcpp::List pseudovarPartial_C(arma::mat Z, Rcpp::List u_indices);
Rcpp::List pseudovarPartial_G(arma::mat Z, const arma::mat& G, Rcpp::List u_indices);
Rcpp::List pseudovarPartial_P(Rcpp::List V_partial, const arma::mat& P);
Rcpp::List pseudovarPartial_V(Rcpp::List V_partial, const arma::mat& V_star_inv);

#endif
