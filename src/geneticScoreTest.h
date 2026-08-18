#ifndef GENETICSCORETEST_H
#define GENETICSCORETEST_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat buildVstar(const arma::vec& wdiag, const arma::mat& Z, const arma::mat& K,
                     const arma::vec& sigmas, const Rcpp::List& u_indices);
Rcpp::List buildVpartial(const arma::mat& Z, const arma::mat& K, const Rcpp::List& u_indices,
                         const int& c);
#endif
