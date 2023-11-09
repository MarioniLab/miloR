#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "multiP.h"
using namespace Rcpp;

List multiP(List partials, arma::mat psvar_in){
    int nsize = partials.size();
    int ps_nrows = psvar_in.n_rows; // this should be the column dimension
    List out(nsize);

    for(int k = 0; k < nsize; k++){
        arma::mat _p = partials[k];
        int ncol = _p.n_cols;

        arma::mat _P(ps_nrows, ncol); // this is an empty sparse matrix rows from psvar_in, cols from _p
        _P = (psvar_in * _p);

        out[k] = _P;
    }

    return out;
}
