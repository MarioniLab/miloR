#include<RcppArmadillo.h>
#include<Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"

// utility functions
// [[Rcpp::export]]
Rcpp::LogicalVector check_na_arma_numeric(arma::vec X){
    // don't being function names with '_'
    // input is an arma::vec
    const int& n = X.size();
    Rcpp::LogicalVector _out(n);
    bool _isnan;

    for(int i=0; i < n; i++){
        // check for NA of any time
        _isnan = std::isnan(X[i]);
        _out[i] = _isnan;
    }

    return _out;
}


// [[Rcpp::export]]
Rcpp::LogicalVector check_inf_arma_numeric(arma::vec X){
    // don't being function names with '_'
    // input is an arma::vec
    const int& n = X.size();
    Rcpp::LogicalVector _out(n);
    bool _isinf;

    for(int i=0; i < n; i++){
        // check for NA of any time
        _isinf = !std::isfinite(X[i]);
        _out[i] = _isinf;
    }

    return _out;
}


// [[Rcpp::export]]
Rcpp::LogicalVector check_zero_arma_numeric(arma::vec X){
    // don't being function names with '_'
    // input is an arma::vec
    const int& n = X.size();
    Rcpp::LogicalVector _out(n);
    bool _iszero;

    for(int i=0; i < n; i++){
        // check for NA of any time
        _iszero = X[i] == 0;
        _out[i] = _iszero;
    }

    return _out;
}



// [[Rcpp::export]]
Rcpp::LogicalVector check_zero_arma_complex(arma::cx_vec X){
    // don't being function names with '_'
    // input is an arma::vec
    const int& n = X.size();
    Rcpp::LogicalVector _out(n);
    bool _iszero;

    for(int i=0; i < n; i++){
        // check for NA of any time
        _iszero = X[i] == 0.0;
        _out[i] = _iszero;
    }

    return _out;
}
