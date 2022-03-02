#ifndef INVERTPSEUDOVAR_H
#define INVERTPSEUDOVAR_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat invertPseudoVar(arma::mat A, arma::mat B, arma::mat Z);

#endif
