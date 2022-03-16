#include "inference.h"
#include "utils.h"
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
    arma::mat _seInv(_se.i());
    arma::vec se(arma::sqrt(_seInv.diag()));

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
    // compute the variance-covariance matrix of the sigma estimates
    //     V_a <- Matrix(0L, nrow=length(curr_sigma), ncol=length(curr_sigma))
    //
    //     for(i in seq_along(V_partial)){
    //         for(j in seq_along(V_partial)){
    // ## broadcast out to the individual RE levels
    //             V_a[i, j] <- 2*(1/(matrix.trace(V_partial[[i]] %*% V_partial[[j]]))) # V_partial contains P * dV/dSigma
    //         }
    //     }

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

///////////////////////////////////////////////////////////////////////////////////////////////
// ---- come back to this as I first need to implement a numerical derivative function ---- //
///////////////////////////////////////////////////////////////////////////////////////////////
// arma::vec computeDF(const arma::mat& X, const arma::vec& se, const arma::mat coeff_mat,
//                     const arma::vec& sigmas, const arma::vec& betas, const arma::mat& P,
//                     Rcpp::List Vpartials){
//     // ###---- first calculate g = derivative of C with respect to sigma ----
//     //     function_jac <- function(x, X.fun=as.matrix(X), W_inv.fun=as.matrix(W_inv), full.Z.fun=as.matrix(full.Z)) {
//     //         UpperLeft <- t(X.fun) %*% W_inv.fun %*% X.fun
//     //         UpperRight <- t(X.fun) %*% W_inv.fun %*% full.Z.fun
//     //         LowerLeft <- t(full.Z.fun) %*% W_inv.fun %*% X.fun
//     //         LowerRight <- t(full.Z.fun) %*% W_inv.fun %*% full.Z.fun
//     //         n <- length(random.levels)
//     //         diag(LowerRight) <- diag(LowerRight) + rep(1/x, times=lengths(random.levels)) #when extending to random slopes, this needs to be changed to a matrix and added to LowerRight directly
//     //         C <- solve(UpperLeft - UpperRight %*% solve(LowerRight) %*% LowerLeft)
//     //     }
//     //
//     //     jac <- jacobian(func=function_jac, x=as.vector(curr_sigma))
//     //     jac_list <- lapply(1:ncol(jac), function(i)
//     //                            array(jac[, i], dim=rep(length(curr_beta), 2))) #when extending to random slopes, this would have to be reformatted into list, where each element belongs to one random effect
//
//     /// This uses the jacobian() function from numDeriv <- which means I'll need to re-implement it in Rcpp
//
//
//     // #next, calculate V_a, the asymptotic covariance matrix of the estimated covariance parameters
//     // #given by formula below
//     //     P <- computeP_REML(V_star_inv=V_star_inv, X=X)
//     //     V_a <- list()
//     //     count <- 1
//     //     for (k in names(random.levels)){
//     //         temp_Vpartial <- V_partial[grep(k, names(V_partial))]
//     //         for (i in 1:length(temp_Vpartial)) {
//     //             for (j in 1:length(temp_Vpartial)) {
//     //                 V_a[[count]] <- 2*(1/(matrix.trace(P %*% temp_Vpartial[[i]] %*% P %*% temp_Vpartial[[j]])))
//     //                 count <- count + 1
//     //             }
//     //         }
//     //     }
//     //     V_a <- bdiag(V_a)
//     //
//     // # P <- computeP_REML(V_star_inv=V_star_inv, X=X)
//     // # V_a <- matrix(NA, nrow=length(random.levels), ncol=length(random.levels))
//     // # for (i in 1:length(random.levels)) {
//     // #   for (j in 1:length(random.levels)) {
//     // #     V_a[i,j] <- 2*1/(matrix.trace(P %*% V_partial[[i]] %*% P %*% V_partial[[j]]))
//     // #   }
//     // # }
//     // # print(V_a)
//     //
//     //         df <- rep(NA, length(curr_beta))
//     //             for (i in 1:length(curr_beta)) {
//     //                 jac_var_beta <- unlist(lapply(lapply(jac_list, diag), `[[`, i))
//     //                 denom <- t(jac_var_beta) %*% (V_a) %*% jac_var_beta #g' Va g
//     //                 df[i] <- 2*((SE[i]^2)^2)/denom
//     //             }
// }
