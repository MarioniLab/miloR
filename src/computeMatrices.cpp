#include "computeMatrices.h"
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec computeYStar(arma::mat X, arma::vec curr_beta, arma::mat Z, arma::mat Dinv, arma::vec curr_u, arma::vec y){
    // compute pseudovariable
    int n = X.n_rows;
    arma::vec ystar(n);

    ystar = ((X * curr_beta) + (Z * curr_u)) + Dinv * (y - exp((X * curr_beta) + (Z * curr_u)));
    return ystar;
}


arma::mat computeVmu(arma::vec mu, double r){
    int n = mu.size();
    arma::mat Vmu(n, n);

    Vmu.diag() = (pow(mu, 2)/r) + mu;

    return Vmu;
}


arma::mat computeW(arma::mat Dinv, arma::mat V){
    int n = Dinv.n_cols;
    arma::mat W(n, n);

    W = Dinv * V * Dinv; // this is a bottle neck - can we speed up these big multiplications?
    return W;
}


arma::mat computeVStar(arma::mat Z, arma::mat G, arma::mat W){
    int n = Z.n_rows;
    arma::mat vstar(n, n);
    vstar = (Z * G * Z.t()) + W;

    return vstar;
}


arma::mat computePREML (arma::mat Vsinv, arma::mat X){
    int n = Vsinv.n_cols;
    // int m = X.n_cols;
    // arma::mat _mid(m , m);
    // _mid = inv(X.t() * Vsinv * X);
    arma::mat P(n, n);
    P = Vsinv - (Vsinv * X * inv(X.t() * Vsinv * X) * X.t() * Vsinv); // also slow with all these multiplications
    return P;
}


arma::mat initialiseG (List u_indices, arma::vec sigmas){
    // construct the correct size of G given the random effects and variance components
    // the independent sigmas go on the diagonal and the off-diagonal are the crossed/interactions
    // this doesn't actually handle the off-diagonal interactions yet
    int c = u_indices.size();
    int stot = 0;

    // sum total number of levels
    for(int i=0; i < c; i++){
        StringVector _ir = u_indices(i);
        stot += _ir.size();
    }

    arma::mat G(stot, stot);
    G = G.zeros();

    // this only fills the diagonal elements of G
    unsigned long i = 0;
    unsigned long j = 0;

    for(unsigned long k=0; k < stot; k++){
        i = k;
        j = k;
        for(int x = 0; x < c; x++){
            arma::uvec _r = u_indices(x);
            unsigned long q = _r.size();
            double _s = sigmas(x);

            for(int l=0; l < q; l++){
                unsigned long _lu = _r(l);

                if(k == _lu - 1){
                    G(i, j) = _s;
                }
            }
        }
    }

    return G;
}


arma::mat invGmat (List u_indices, arma::vec sigmas){
    // first construct the correct sized G, i.e. c x c, then brodcast this to all RE levels
    // make little G inverse
    int c = u_indices.size();
    int stot = 0;
    arma::vec lsigma(c);

    for(int k = 0; k < c; k++){
        lsigma(k) = 1/sigmas(k);
    }

    // sum total number of levels
    for(int i=0; i < c; i++){
        StringVector _ir = u_indices(i);
        stot += _ir.size();
    }

    arma::mat G(stot, stot);
    G = G.zeros();

    // this only fills the diagonal elements of G
    unsigned long i = 0;
    unsigned long j = 0;

    for(unsigned long k=0; k < stot; k++){
        i = k;
        j = k;
        for(int x = 0; x < c; x++){
            arma::uvec _r = u_indices(x);
            unsigned long q = _r.size();
            double _s = lsigma(x);

            for(int l=0; l < q; l++){
                unsigned long _lu = _r(l);

                if(k == _lu - 1){
                    G(i, j) = _s;
                }
            }
        }
    }

    return G;
}




