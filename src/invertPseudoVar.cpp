#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "invertPseudoVar.h"
using namespace Rcpp;

arma::mat invertPseudoVar(arma::mat A, arma::mat B, arma::mat Z){
    // make sparse - helps with the many matrix multiplications
    int c = B.n_cols;
    int n = A.n_cols;

    arma::sp_mat spA(A);
    arma::sp_mat spB(B);
    arma::sp_mat spZ(Z);

    arma::sp_mat I = arma::speye<arma::sp_mat>(c, c); // create the cxc identity matrix
    arma::sp_mat omt(n, n);
    arma::sp_mat mid(c, c);
    arma::sp_mat AZB(A.n_rows, B.n_cols);
    arma::mat out_omt(omt);

    AZB = spA * spZ * spB;
    mid = I + (spZ.t() * AZB); // If we know the structure in B can we simplify this more???
    arma::mat dmid(mid);

    try{
        // double _rcond = arma::rcond(mid);
        double _rcond = arma::rcond(dmid);
        bool is_singular;
        is_singular = _rcond < 1e-12;

        // check for singular condition
        if(is_singular){
            Rcpp::stop("Pseudovariance component matrix is computationally singular");
        }

        arma::sp_mat midinv(dmid.i());
        omt = spA - (AZB * midinv * spZ.t() * spA); // stack multiplications like this appear to be slow
        arma::mat out_omt(omt);
        return out_omt;
    } catch(std::exception &ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }

    return out_omt;
}
