#include<RcppArmadillo.h>
#include<RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include "paramEst.h"
#include "computeMatrices.h"
#include "invertPseudoVar.h"
#include "pseudovarPartial.h"
#include "multiP.h"
#include "inference.h"
using namespace Rcpp;

//' GLMM parameter estimation using pseudo-likelihood
//'
//' Iteratively estimate GLMM fixed and random effect parameters, and variance
//' component parameters using Fisher scoring based on the Pseudo-likelihood
//' approximation to a Normal loglihood.
//' @param Z sp_mat - sparse matrix that maps random effect variable levels to
//' observations
//' @param X sp_mat - sparse matrix that maps fixed effect variables to
//' observations
//' @param muvec NumericVector vector of estimated phenotype means
//' @param curr_theta NumericVector vector of initial parameter estimates
//' @param curr_beta NumericVector vector of initial beta estimates
//' @param curr_u NumericVector of initial u estimates
//' @param curr_G NumericVector of initial sigma estimates
//' @param y NumericVector of observed counts
//' @param rlevels List containing mapping of RE variables to individual
//' levels
//' @param curr_disp double Dispersion parameter estimate
//' @param REML bool - use REML for variance component estimation
// [[Rcpp::export]]
List fitPLGlmm(const arma::mat& Z, const arma::mat& X, arma::vec muvec, arma::vec curr_beta,
               arma::vec curr_theta, arma::vec curr_u, arma::vec curr_sigma,
               arma::mat curr_G, const arma::vec& y, List u_indices,
               double theta_conv,
               const List& rlevels, double curr_disp, const bool& REML, const int& maxit){

    // declare all variables
    List outlist(10);
    int iters=0;
    int stot = Z.n_cols;
    const int& c = curr_sigma.size();
    const int& m = X.n_cols;
    const int& n = X.n_rows;
    bool meet_cond = false;

    // setup matrices
    arma::mat D(n, n);
    D.zeros();
    arma::mat Dinv(n, n);
    Dinv.zeros();

    arma::vec y_star(n);

    arma::mat Vmu(n, n);
    arma::mat W(n, n);
    arma::mat Winv(n, n);

    arma::mat V_star(n, n);
    arma::mat V_star_inv(n, n);
    arma::mat P(n, n);

    arma::mat coeff_mat(m+c, m+c);
    List V_partial(c);

    arma::vec score_sigma(c);
    arma::mat information_sigma(c, c);
    arma::vec sigma_update(c);
    arma::vec sigma_diff(sigma_update.size());
    sigma_diff.zeros();

    arma::mat G_inv(stot, stot);

    arma::vec theta_update(m+stot);
    arma::vec theta_diff(theta_update.size());
    theta_diff.zeros();

    // setup vectors to index the theta updates
    // assume always in order of beta then u
    arma::uvec beta_ix(m);
    for(int x=0; x < m; x++){
        beta_ix[x] = x;
    }

    arma::uvec u_ix(stot);
    for(int px = 0; px < stot; px++){
        u_ix[px] = m + px;
    }

    bool converged = false;
    while(!meet_cond){
        D.diag() = muvec;
        Dinv = D.i();
        y_star = computeYStar(X, curr_beta, Z, Dinv, curr_u, y);
        Vmu = computeVmu(muvec, curr_disp);
        W = computeW(Dinv, Vmu);
        Winv = W.i();
        V_star = computeVStar(Z, curr_G, W);
        V_star_inv = invertPseudoVar(Winv, curr_G, Z);

        if(REML){
            P = computePREML(V_star_inv, X);
            V_partial = pseudovarPartial_P(Z, u_indices, P);

            score_sigma = sigmaScoreREML_arma(V_partial, y_star, P);
            information_sigma = sigmaInfoREML_arma(V_partial, P);
        } else{
            V_partial = pseudovarPartial_C(Z, u_indices);
            score_sigma = sigmaScore(y_star, curr_beta, X, V_partial, V_star_inv);
            information_sigma = sigmaInformation(V_star_inv, V_partial);
        };

        sigma_update = FisherScore(information_sigma, score_sigma, curr_sigma);
        sigma_diff = abs(sigma_update - curr_sigma); // needs to be an unsigned real value

        // update sigma, G, and G_inv
        curr_sigma = sigma_update;
        curr_G = initialiseG(u_indices, curr_sigma);
        G_inv = invGmat(u_indices, curr_sigma);

        // Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u
        // compute the coefficient matrix
        coeff_mat = coeffMatrix(X, Winv, Z, G_inv);
        theta_update = solve_equations(stot, m, Winv, Z.t(), X.t(), coeff_mat, curr_beta, curr_u, y_star);
        theta_diff = abs(theta_update - curr_theta);

        curr_theta = theta_update;
        curr_beta = curr_theta.elem(beta_ix); // how do we subset the correct elements without having names? Need a record of the indicies
        curr_u = curr_theta.elem(u_ix);

        muvec = exp((X * curr_beta) + (Z * curr_u));
        iters++;

        bool _thconv = false;
        _thconv = all(theta_diff < theta_conv);

        bool _siconv = false;
        _siconv = all(sigma_diff < theta_conv);

        bool _ithit = false;
        _ithit = iters > maxit;

        meet_cond = ((_thconv && _siconv) || _ithit);
        converged = _thconv && _siconv;
    }

    arma::vec se(computeSE(m, c, coeff_mat));
    arma::vec tscores(computeTScore(curr_beta, se));
    arma::mat vcov(varCovar(V_partial, c));
    // It will be expensive to sequentially grow a list here - do we even need to return
    // all of the intermediate results??
    outlist = List::create(_["FE"]=curr_beta, _["RE"]=curr_u, _["Sigma"]=curr_sigma,
                           _["converged"]=converged, _["Iters"]=iters, _["Dispersion"]=curr_disp,
                           _["Hessian"]=information_sigma, _["SE"]=se, _["t"]=tscores,
                           _["COEFF"]=coeff_mat, _["P"]=P, _["Vpartial"]=V_partial, _["Ginv"]=G_inv,
                           _["Vsinv"]=V_star_inv, _["Winv"]=Winv, _["VCOV"]=vcov);

    return outlist;
}


