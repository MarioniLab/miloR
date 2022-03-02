#include "computeMatrices.h"
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec computeYStar(arma::mat X, arma::vec curr_beta, arma::mat Z, arma::mat Dinv, arma::vec curr_u, arma::vec y){
    // compute pseudovariable
    int n = X.n_rows;
    arma::vec ystar(n);

    ystar = ((X * curr_beta) + (Z * curr_u)) + Dinv * (y - exp((X * curr_beta) + (Z * curr_u)));
    return ystar;
}


// [[Rcpp::export]]
arma::mat computeVmu(arma::vec mu, double r){
    int n = mu.size();
    arma::mat Vmu(n, n);

    Vmu.diag() = (pow(mu, 2)/r) + mu;

    return Vmu;
}


// [[Rcpp::export]]
arma::mat computeW(arma::mat Dinv, arma::mat V){
    int n = Dinv.n_cols;
    arma::mat W(n, n);

    W = Dinv * V * Dinv;
    return W;
}


// [[Rcpp::export]]
arma::mat computeVStar(arma::mat Z, arma::mat G, arma::mat W){
    int n = Z.n_rows;
    arma::mat vstar(n, n);
    vstar = (Z * (G * Z.t())) + W;

    return vstar;
}


// [[Rcpp::export]]
arma::mat computePREML (arma::mat Vsinv, arma::mat X){
    // Vs_inv -  Vsinv * X * [t(X) * Vs_inv * X]^-1 * t(X) * Vsinv
    int n = Vsinv.n_cols;
    arma::mat P(n, n);

    P = Vsinv - (Vsinv * X * inv(X.t() * Vsinv * X) * X.t() * Vsinv);
    return P;
}


// [[Rcpp:export]]
arma::mat initialiseG (List rlevels, arma::vec sigmas){
    // construct the correct size of G given the random effects and variance components
    // the independent sigmas go on the diagonal and the off-diagonal are the crossed/interactions
    // this doesn't actually handle the off-diagonal interactions yet

    int c = rlevels.size();
    int stot = 0;

    // sum total number of levels
    for(int i=0; i < c; i++){
        arma::vec _ir = rlevels[i];
        stot += _ir.size();
    }

    arma::mat G(stot, stot);
    G = G.zeros();

    int i = 0;
    int j = 0;
    for(int k=0; k < c; k++){
        arma::vec _r = rlevels[k];
        int q = _r.size();
        double _s = sigmas[k];

        for(int l=0; l < q; l++){
            for(int x=0; x < q; x++){
                G[i, j] = _s;
                i += l;
                j += x;
            }
        }
    }

    return G;
}



