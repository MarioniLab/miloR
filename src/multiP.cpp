#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "multiP.h"
using namespace Rcpp;

//' Compute product of each pseudovariance partial derivatives with the inverse
//' pseudovariance matrix
//'
//' For each variance component, we compute the matrix multiplication of the
//' relevant partial derivative of dV_start/dSigm with the pseudovariance matrix
//'
//' @param List partials - list containing matrices of partial derivatives of the pseudovariance
//' for each variance component
//' @param SparseMatrix psvar_in - inverse of the pseudovariance matrix
// [[Rcpp::export]]
List multiP(List partials, arma::mat psvar_in){
    // each list component is a sparse matrix
    // psvar_in cannot be a dgeMatrix
    // this isn't any faster than the equivalent Rcode - why?
    // Is all the interconversion creating a new bottleneck? I don't think so.
    int nsize = partials.size();
    int ps_nrows = psvar_in.n_rows; // this should be the column dimension
    List out(nsize);

    for(int k = 0; k < nsize; k++){
        arma::sp_mat _p = partials[k];
        int ncol = _p.n_cols;

        arma::sp_mat _P(ps_nrows, ncol); // this is an empty sparse matrix rows from psvar_in, cols from _p
        _P = (psvar_in * _p);

        out[k] = _P;
    }

    return out;
}
