#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
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
    int nsize = partials.size();
    int ps_nrows = psvar_in.n_rows; // this should be the column dimension
    List out(nsize);

    for(int k = 0; k < nsize; k++){
        S4 _p = partials[k];

        // get the sparse matrix elements from the S4 object to make and arma::mat
        IntegerVector dims = _p.slot("Dim");
        arma::urowvec i = Rcpp::as<arma::urowvec>(_p.slot("i"));
        arma::urowvec p = Rcpp::as<arma::urowvec>(_p.slot("p"));
        arma::vec x = Rcpp::as<arma::vec>(_p.slot("x"));

        int nrow = dims[0], ncol = dims[1];
        arma::sp_mat imat(i, p, x, nrow, ncol);

        arma::mat _P(ps_nrows, ncol); // this is an empty sparse matrix rows from psvar_in, cols from _p
        _P = (psvar_in * imat);

        out[k] = _P;
    }

    return out;
}
