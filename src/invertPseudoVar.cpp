#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Compute the inverse of a structured covariance matrix
//'
//' Using Henderson's adjusted Woodbury formula for a singular B matrix,
//' compute the inverse of the pseudocovariance matrix ZGZ' + W as
//' (A + UBU^T)^-1 = A^-1 - A^-1UB[I + U^TA^-1UB]^-1U^TA^-1
//'
//' @param A SparseMatrix - a nxn matrix of the GLMM covariance D^-1*V*D^-1
//' @param B SparseMatrix - a cxc matrix of variance components
//' @param Z SparseMatrix - a nxc design matrix that maps REs to samples
// [[Rcpp::export]]
arma::sp_mat invertPseudoVar(arma::sp_mat A, arma::sp_mat B, arma::sp_mat Z){
    int c = B.n_cols;
    int n = A.n_cols;
    arma::sp_mat I = arma::eye<arma::sp_mat>(c, c); // create the cxc identity matrix
    arma::sp_mat omt(n, n);
    arma::mat mid(c, c);
    mid = I + Z.t() * A * Z * B;
    arma::mat midinv(c, c);
    midinv = mid.i();

    omt = A * Z * B * midinv * Z.t() * A;
    return omt;
}
