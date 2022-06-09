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


arma::mat computeW(double disp, arma::mat Dinv, arma::mat V){
    int n = Dinv.n_cols;
    arma::mat W(n, n);

    // this should be trivial as it is a diagonal matrix
    // D^-1 * V_mu * D^-1 simplifies to a diagonal matrix
    // of disp + 1/mu_i which is (phi * I) + Dinv <- we don't need any multiplication!!
    arma::fmat idisp = arma::eye<arma::fmat>(n, n);
    idisp = (1/disp) * idisp;
    W = idisp + Dinv;
    // W = Dinv * V * Dinv; // this is a bottle neck - can we speed up these big multiplications?
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


arma::mat subMatG (double sigma, arma::mat broadcast){
    // construct the submatrix for the input variance component
    int nrow = broadcast.n_rows;
    int ncol = broadcast.n_cols;

    arma::mat subG(nrow, ncol);
    subG = sigma * broadcast;

    return subG;
}



arma::mat initialiseG_G (List u_indices, arma::vec sigmas, arma::mat Kin){
    // construct the correct size of G given the random effects and variance components
    // the independent sigmas go on the diagonal and the off-diagonal are the crossed/interactions
    // this doesn't actually handle the off-diagonal interactions yet
    // for the arbitrary covariance case multiply by Kin
    // is the "genetic" sigma always last?
    int c = u_indices.size();
    int stot = 0;

    // sum total number of levels
    for(int i=0; i < c; i++){
        StringVector _ir = u_indices(i);
        stot += _ir.size();
    }

    // I need to construct the full G from all of the sub-matrices
    // there will be c submatrices to create block-diagonal matrix
    // with zeros on the off-diagonals

    // construct the big matrix for each pair
    // start with the first matrix - then grow a new one at each iteration
    // arma::mat G(stot, stot);
    // G = G.zeros();

    // // this only fills the diagonal elements of G
    // unsigned long i = 0;
    // unsigned long j = 0;
    //
    // for(unsigned long k=0; k < stot; k++){
    //     i = k;
    //     j = k;
    //     for(int x = 0; x < c; x++){
    //         arma::uvec _r = u_indices(x); // the vector of indices of Z that map to the RE
    //         unsigned long q = _r.size(); // the number of levels for the RE
    //         double _s = sigmas(x); // the sigma of the RE
    //
    //         // flow control, if the last effect then it must be the kin ship matrix
    //         if(x == c - 1){
    //             unsigned long n = Kin.n_cols;
    //
    //             if(q != n){
    //                 stop("RE indices and dimensions of covariance do not match");
    //             } else{
    //                 double _sg = sigmas(x);
    //
    //                 for(int l=0; l < n; l++){
    //                     unsigned long _lu = _r(l);
    //
    //
    //                     if(k == _lu - 1){
    //                         G(i, j) = _sg;
    //                     }
    //                 }
    //             }
    //         } else{
    //             for(int l=0; l < q; l++){
    //                 unsigned long _lu = _r(l);
    //                 if(k == _lu - 1){
    //                     G(i, j) = _s;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // return G;
    arma::uvec _Gindex(c); // G is always square
    Rcpp::List Glist(1); // to store the first G

    for(int x = 0; x < c; x++){
        // create the broadcast matrix
        arma::uvec _r = u_indices(x); // the vector of indices of Z that map to the RE
        unsigned long q = _r.size(); // the number of levels for the RE
        double _s = sigmas(x); // the sigma of the RE

        arma::mat sG(q, q);

        if(x == c - 1){
            unsigned long n = Kin.n_cols;
            if(q != n){
                stop("RE indices and dimensions of covariance do not match");
            } else{
                sG = subMatG(_s, Kin);
            }
        } else{
            // create the rxr identity matrix
            arma::mat rEye(q, q, arma::fill::eye);
            sG = subMatG(_s, rEye);
        }

        // grow G at each iteration here
        if(x == 0){
            unsigned long ig_cols = sG.n_cols;
            unsigned long ig_rows = sG.n_rows;
            Glist(0) = sG;
        } else{
            unsigned long sg_cols = sG.n_cols;
            unsigned long sg_rows = sG.n_rows;

            arma::mat G = Glist(0);

            unsigned long g_cols = G.n_cols;
            unsigned long g_rows = G.n_rows;

            arma::mat gright(g_rows, sg_cols);
            arma::mat gleft(sg_rows, g_cols);

            arma::mat top(g_rows, g_cols + sg_cols);
            arma::mat bottom(sg_rows, sg_cols + g_cols);

            top = arma::join_rows(G, gright);
            bottom = arma::join_rows(gleft, sG);

            arma::mat _G(sg_rows + g_rows, sg_cols + g_cols);
            _G = arma::join_cols(top, bottom);
            Glist(0) = _G;
        }
    }

    arma::mat G = Glist(0);
    return G;
}


arma::mat invGmat_G (List u_indices, arma::vec sigmas, arma::mat Kin){
    // first construct the correct sized G, i.e. c x c, then brodcast this to all RE levels
    // make little G inverse
    int c = u_indices.size();
    int stot = 0;

    // sum total number of levels
    for(int i=0; i < c; i++){
        StringVector _ir = u_indices(i);
        stot += _ir.size();
    }

    // arma::mat invG(stot, stot);
    // invG = G.zeros();

    arma::uvec _Gindex(c); // G is always square
    Rcpp::List Glist(1); // to store the first G

    for(int x = 0; x < c; x++){
        // create the broadcast matrix
        arma::uvec _r = u_indices(x); // the vector of indices of Z that map to the RE
        unsigned long q = _r.size(); // the number of levels for the RE
        double _s = sigmas(x); // the sigma of the RE
        double _sinv = 1/_s;

        arma::mat sG(q, q);

        if(x == c - 1){
            unsigned long n = Kin.n_cols;
            if(q != n){
                stop("RE indices and dimensions of covariance do not match");
            } else{
                sG = subMatG(_sinv, Kin); // this needs to be Kin^-1
            }
        } else{
            // create the rxr identity matrix
            arma::mat rEye(q, q, arma::fill::eye);
            sG = subMatG(_sinv, rEye);
        }

        // grow G at each iteration here
        if(x == 0){
            unsigned long ig_cols = sG.n_cols;
            unsigned long ig_rows = sG.n_rows;
            Glist(0) = sG;
        } else{
            unsigned long sg_cols = sG.n_cols;
            unsigned long sg_rows = sG.n_rows;

            arma::mat G = Glist(0);

            unsigned long g_cols = G.n_cols;
            unsigned long g_rows = G.n_rows;

            arma::mat gright(g_rows, sg_cols);
            arma::mat gleft(sg_rows, g_cols);

            arma::mat top(g_rows, g_cols + sg_cols);
            arma::mat bottom(sg_rows, sg_cols + g_cols);

            top = arma::join_rows(G, gright);
            bottom = arma::join_rows(gleft, sG);

            arma::mat _G(sg_rows + g_rows, sg_cols + g_cols);
            _G = arma::join_cols(top, bottom);
            Glist(0) = _G;
        }
    }

    arma::mat G = Glist(0);

    // // this only fills the diagonal elements of G
    // unsigned long i = 0;
    // unsigned long j = 0;
    //
    // for(unsigned long k=0; k < stot; k++){
    //     i = k;
    //     j = k;
    //     for(int x = 0; x < c; x++){
    //         arma::uvec _r = u_indices(x);
    //         unsigned long q = _r.size();
    //         double _s = lsigma(x);
    //
    //         for(int l=0; l < q; l++){
    //             unsigned long _lu = _r(l);
    //
    //             if(k == _lu - 1){
    //                 G(i, j) = _s;
    //             }
    //         }
    //     }
    // }
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




