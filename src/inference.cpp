#include "inference.h"
#include<RcppArmadillo.h>
#include<RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// using namespace Rcpp;

// All functions used for inference

arma::vec computeSE(const int& m, const int& c, const arma::mat& coeff_mat) {
    // compute the fixed effect standard errors from the MME coefficient matrix
    arma::mat ul(m, m);
    ul = coeff_mat(arma::span(0, m-1), arma::span(0, m-1)); // m X m
    arma::mat ur(coeff_mat.submat(0, m, m-1, m+c-1)); // m X c
    arma::mat ll(coeff_mat.submat(m, 0, m+c-1, m-1)); // c X m
    arma::mat lr(coeff_mat.submat(m, m, m+c-1, m+c-1)); // c X c

    arma::mat _se(ul - ur * lr.i() * ll); // m X m - (m X c X m) <- this should commute
    // will need a check here for singular hessians...
    try{
        double _rcond = arma::rcond(_se);
        bool is_singular;
        is_singular = _rcond < 1e-12;

        // check for singular condition
        if(is_singular){
            Rcpp::stop("Standard Error coefficient matrix is computationally singular");
        }

        arma::mat _seInv(_se.i());
        arma::vec se(arma::sqrt(_seInv.diag()));
        return se;
    } catch(std::exception &ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }


}


arma::vec computeTScore(const arma::vec& curr_beta, const arma::vec& SE){

    const int& m = curr_beta.size();
    const int& selength = SE.size();

    if(m != selength){
        Rcpp::stop("standard errors and beta estimate sizes differ");
    }

    arma::vec tscore(m);

    for(int x=0; x < m; x++){
        double _beta = curr_beta[x];
        double _se = SE[x];

        double _t = _beta/_se;
        tscore[x] = _t;
    }

    return tscore;
}


arma::mat varCovar(const Rcpp::List& psvari, const int& c){

    arma::mat Va(c, c);
    for(int i=0; i < c; i++){
        arma::mat _ips = psvari(i);
        for(int j=i; j < c; j++){
            arma::mat _jps = psvari(j);
            arma::mat _ij(_ips * _jps);
            Va(i, j) = 2 * (1/(arma::trace(_ij)));
            if(i != j){
                Va(j, i) = 2 * (1/(arma::trace(_ij)));
            }
        }
    }

    return Va;
}

