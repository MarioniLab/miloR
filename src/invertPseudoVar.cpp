#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "invertPseudoVar.h"
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
arma::mat invertPseudoVar(arma::mat A, arma::mat B, arma::mat Z){
    int c = B.n_cols;
    int n = A.n_cols;
    arma::mat I = arma::eye<arma::mat>(c, c); // create the cxc identity matrix
    arma::mat omt(n, n);
    arma::mat mid(c, c);
    arma::mat AZB(A.n_rows, B.n_cols);
    AZB = A * Z * B;
    mid = I + (Z.t() * AZB); // If we know the structure in B can we simplify this more???
    arma::mat midinv(c, c);

    try{
        double _rcond = arma::rcond(mid);
        bool is_singular;
        is_singular = _rcond < 1e-12;

        // check for singular condition
        if(is_singular){
            Rcpp::stop("Pseudovariance component matrix is computationally singular");
        }

        midinv = mid.i();
        omt = A - (AZB * midinv * Z.t() * A); // stack multiplications like this appear to be slow
        return omt;
    } catch(std::exception &ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }
}
