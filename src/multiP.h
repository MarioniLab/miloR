#ifndef MULTIP_H
#define MULTIP_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List multiP(Rcpp::List partials, arma::mat psvar_in);

#endif
