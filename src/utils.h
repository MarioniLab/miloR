#ifndef UTILS_H
#define UTILS_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::LogicalVector check_na_arma_numeric(arma::vec x);
Rcpp::LogicalVector check_inf_arma_numeric(arma::vec X);
Rcpp::LogicalVector check_zero_arma_numeric(arma::vec X);
#endif
