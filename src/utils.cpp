#include<RcppArmadillo.h>
#include<Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"

// utility functions
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


Rcpp::LogicalVector check_zero_arma_complex(arma::cx_vec X){
    // don't begin function names with '_'
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


Rcpp::LogicalVector check_tol_arma_numeric(arma::vec X, double tol){
    // don't being function names with '_'
    // input is an arma::vec
    const int& n = X.size();
    Rcpp::LogicalVector _out(n);
    bool _iszero;

    for(int i=0; i < n; i++){
        // check for NA of any time
        _iszero = X[i] <= tol;
        _out[i] = _iszero;
    }

    return _out;
}

bool check_pd_matrix(arma::mat A){
    // check that A matrix is positive definite - i.e. all positive eigenvalues
    // A must be square
    unsigned int m = A.n_cols;
    unsigned int n = A.n_rows;
    bool _is_sym;

    if(m != n){
        return false;
    }

    _is_sym = A.is_symmetric();

    if(!_is_sym){
        Rcpp::stop("matrix A is not symmetric");
        return false;
    }

    arma::vec eigenvals = arma::eig_sym(A);
    bool alltrue;
    alltrue = all(eigenvals > 0.0);
    return alltrue;
}
