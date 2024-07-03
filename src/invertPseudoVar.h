#ifndef INVERTPSEUDOVAR_H
#define INVERTPSEUDOVAR_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat invertPseudoVar(const arma::mat& A, const arma::mat& B, const arma::mat& Z,
                          const arma::mat& ZtA);
arma::mat kRankOneUpdates(const arma::mat& Vinv, const arma::mat& B);
arma::mat rankOneUp(const arma::mat& A, const arma::uvec& u, const arma::drowvec& v);
#endif
