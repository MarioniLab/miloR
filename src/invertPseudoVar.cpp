#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "invertPseudoVar.h"
using namespace Rcpp;

arma::mat invertPseudoVar(arma::mat A, arma::mat B, arma::mat Z){
    int c = B.n_cols;
    int n = A.n_cols;
    arma::mat I = arma::eye<arma::mat>(c, c); // create the cxc identity matrix
    arma::mat omt(n, n, arma::fill::zeros);
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

    return omt;
}
