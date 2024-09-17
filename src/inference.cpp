#include "inference.h"
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;

// All functions used for inference
arma::vec computeSE(const int& m, const int& c, const arma::mat& coeff_mat) {
    // compute the fixed effect standard errors from the MME coefficient matrix

    const int& l = coeff_mat.n_cols; // this should be m + c
    const int& p = coeff_mat.n_rows; // this should be m + c

    const int& scol = m + c;
    const int& srow = m + c;
    if(l != scol){
        Rcpp::Rcout << scol << std::endl;
        Rcpp::Rcout << l << std::endl;
        Rcpp::Rcout << p << std::endl;
        Rcpp::stop("N cols and input dimensions m + c are not equal");
    }

    if(p != srow){
        Rcpp::stop("N rows and input dimensions m + c are not equal");
    }

    arma::mat ul(coeff_mat.submat(0, 0, m-1, m-1)); // m X m
    arma::mat ur(coeff_mat.submat(0, m, m-1, m+c-1)); // m X l
    arma::mat ll(coeff_mat.submat(m, 0, m+c-1, m-1)); // p X m
    arma::mat lr(coeff_mat.submat(m, m, m+c-1, m+c-1)); // p X l


    arma::mat _se(ul - ur * lr.i() * ll); // m X m - (m X c X m) <- this should commute
    arma::vec se(m+c);
    // will need a check here for singular hessians...
    double _rcond = arma::rcond(_se);
    if(_rcond < 1e-12){
        Rcpp::warning("Standard Error coefficient matrix is computationally singular - using pseudoinverse");
        arma::mat _seInv = arma::pinv(_se);
        se = arma::sqrt(_seInv.diag());
    } else{
        arma::mat _seInv = arma::inv(_se);
        se = arma::sqrt(_seInv.diag());
    }

    return se;
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
    // trace can be computed as sum of dot products of rows and columns

    for(int i=0; i < c; i++){
        arma::mat _ips = psvari[i]; // why isn't this
        for(int j=i; j < c; j++){
            arma::mat _jps = psvari(j);

            // Compute trace directly
            double trace = 0.0;
            for(arma::uword k = 0; k < _ips.n_rows; ++k) {
                trace += arma::dot(_ips.row(k), _jps.col(k));
            }

            Va(i, j) = 2 * (1/trace);
            if(i != j){
                Va(j, i) = 2 * (1/trace);
            }
        }
    }

    return Va;
}

