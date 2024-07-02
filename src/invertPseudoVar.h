#ifndef INVERTPSEUDOVAR_H
#define INVERTPSEUDOVAR_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat invertPseudoVar(arma::mat A, arma::mat B, arma::mat Z);
arma::mat kRankOneUpdates(const arma::mat& Vinv, const arma::mat& B);
arma::mat rankOneUp(const arma::mat& A, const arma::uvec& u, const arma::drowvec& v);
#endif
