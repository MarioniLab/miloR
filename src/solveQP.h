#ifndef SOLVEQP_H
#define SOLVEQP_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec solveQP(arma::mat A_arma, arma::vec b_arma, arma::vec x_arma);
#endif
