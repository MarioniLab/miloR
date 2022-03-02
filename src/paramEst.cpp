#include "paramEst.h"
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// All functions used in parameter estimation

// [[Rcpp::export]]
arma::vec sigmaScoreREML (List pvstar_i, arma::mat Vsinv, arma::vec ystar, arma::mat P){

    int c = pvstar_i.size();
    arma::vec reml_score(c);

    for(int i=0; i < c; i++){
        arma::mat _pvi = pvstar_i[i];
        double lhs = -0.5 * arma::trace(_pvi);
        // double mid1 = 0.0;
        arma::mat mid1(1, 1);
        mid1 = arma::trans(ystar) * _pvi * P * ystar;
        mid1.print();
        double rhs = 0.5 * mid1(0, 0);

        reml_score[i] = lhs + rhs;
    }

    return reml_score;
}


// [[Rcpp::export]]
arma::mat sigmaInfoREML (List pvstari){
    // REML Fisher/expected information matrix

    int c = pvstari.size();
    arma::mat sinfo(c, c);

    for(int i=0; i < c; i++){
        arma::mat _ip = pvstari[i];
        for(int j=0; j < c; j++){
            arma::mat _jp = pvstari[j];
            arma::mat _ij(1, 1);
            _ij = _ip * _jp;
            sinfo[i, j] = _ij[0, 0];
        }
    }

    return sinfo;
}


// [[Rcpp::export]]
arma::vec sigmaScore (arma::vec ystar, arma::vec beta, arma::mat X, List V_partial, arma::mat V_star_inv){

    int c = V_partial.size();
    int n = X.n_rows;
    arma::vec score(c);
    arma::vec ystarminx(n);
    ystarminx = ystar - (X * beta);


    for(int i=0; i < c; i++){
        arma::mat _ip = V_partial[i];
        double lhs = -0.5 * arma::trace(V_star_inv * _ip);
        arma::mat rhs_mat(1, 1);

        rhs_mat = ystarminx.t() * V_star_inv * _ip * V_star_inv * ystarminx;
        score[i] = lhs + 0.5 * rhs_mat[0, 0];
    }

    return score;
}


// [[Rcpp::export]]
arma::mat sigmaInformation (arma::mat V_star_inv, List V_partial){
    int c = V_partial.size();
    int n = V_star_inv.n_cols;
    arma::mat sinfo = arma::zeros(c, c);

    for(int i=0; i < c; i++){
        arma::mat _ip = V_partial[i];
        for(int j=0; j < c; j++){
            arma::mat _jp = V_partial[j];

            arma::mat _inmat(n, n);
            double _tr = 0.0;
            _inmat = V_star_inv * _ip * V_star_inv * _jp;
            _tr = 0.5 * arma::trace(_inmat);

            sinfo[i, j] = _tr;
        }
    }

    return sinfo;
}


// [[Rcpp::export]]
arma::vec FisherScore (arma::mat hess, arma::vec score_vec, arma::vec theta_hat){
    // sequentially update the parameter using the Newton-Raphson algorithm
    // theta ~= theta_hat + hess^-1 * score
    // this needs to be in a direction of descent towards a minimum

    int m = theta_hat.size();
    arma::vec theta(m);

    theta = theta_hat + (hess.i() * score_vec);
    return theta;
}


// [[Rcpp::export]]
arma::vec solve_equations (arma::mat X, arma::mat Winv, arma::mat Z, arma::mat Ginv, arma::vec beta, arma::vec u, arma::vec ystar){
    // solve the mixed model equations

    int c = Z.n_cols;
    int m = X.n_cols;

    arma::mat ul(m, m);
    arma::mat ur(m, c);
    arma::mat ll(c, m);
    arma::mat lr(c, c);

    arma::mat lhs_top(m, m+c);
    arma::mat lhs_bot(c, m+c);
    arma::mat lhs(m+c, m+c);

    arma::vec rhs_beta(m);
    arma::vec rhs_u(c);
    arma::vec rhs(m+c);

    arma::vec theta_up(m+c);

    ul = X.t() * Winv * X;
    ur = X.t() * Winv * Z;
    ll = Z.t() * Winv * X;
    lr = (Z.t() * Winv * Z) + Ginv;

    lhs_top = arma::join_cols(ul, ur);
    lhs_bot = arma::join_cols(ll, lr);
    lhs = arma::join_rows(lhs_top, lhs_bot);

    rhs_beta = X.t() * Winv * ystar;
    rhs_u = Z.t() * Winv * ystar;
    rhs = arma::join_rows(rhs_beta, rhs_u);

    theta_up = arma::inv(lhs) * rhs;
    return theta_up;
}

