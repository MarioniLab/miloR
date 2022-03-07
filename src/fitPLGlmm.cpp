#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "paramEst.h"
#include "computeMatrices.h"
#include "invertPseudoVar.h"
#include "pseudovarPartial.h"
#include "multiP.h"
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
List fitPLGlmm(arma::mat Z, arma::mat X, arma::vec muvec, arma::vec curr_beta,
               arma::vec curr_theta, arma::vec curr_u, arma::vec curr_sigma,
               arma::mat curr_G, arma::vec y, List u_indices,
               double theta_conv,
               List rlevels, double curr_disp, bool REML, int maxit){

    // declare all variables
    List outlist(7);
    int iters=0;
    int stot = Z.n_cols;
    int c = curr_sigma.size();
    int m = X.n_cols;
    int n = X.n_rows;
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
        V_partial = pseudovarPartial_C(Z, u_indices);

        if(REML){
            arma::mat P(n, n);
            P = computePREML(V_star_inv, X);
            List PVSTARi(c);
            PVSTARi = multiP(V_partial, V_star_inv);

            score_sigma = sigmaScoreREML(PVSTARi, V_star_inv, y_star, P);
            information_sigma = sigmaInfoREML(PVSTARi);
        } else{
            score_sigma = sigmaScore(y_star, curr_beta, X, V_partial, V_star_inv);
            information_sigma = sigmaInformation(V_star_inv, V_partial);
        };
        score_sigma.print("Sigma score vector");
        information_sigma.print("Fisher Information");

        sigma_update = FisherScore(information_sigma, score_sigma, curr_sigma);
        sigma_diff = abs(sigma_update - curr_sigma); // needs to be an unsigned real value

        // update sigma, G, and G_inv
        curr_sigma = sigma_update;
        curr_G = initialiseG(rlevels, curr_sigma);
        G_inv = curr_G.i(); // is this ever singular?

        // Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u
        theta_update = solve_equations(X, Winv, Z, G_inv, curr_beta, curr_u, y_star);
        theta_diff = abs(theta_update - curr_theta);

        curr_theta = theta_update;
        curr_beta = curr_theta.elem(beta_ix); // how do we subset the correct elements without having names? Need a record of the indicies
        curr_u = curr_theta.elem(u_ix);

        muvec = exp((X * curr_beta) + (Z * curr_u));
        iters++;

        theta_diff.print("theta_diff");
        sigma_diff.print("sigma_diff");
        Rprintf("%0.8f\n", theta_conv);

        bool _thconv = false;
        _thconv = all(theta_diff < theta_conv);
        Rprintf("Theta converge %s\n", _thconv ? "true" : "false");

        bool _siconv = false;
        _siconv = all(sigma_diff < theta_conv);
        Rprintf("Sigma converge %s\n", _siconv ? "true" : "false");

        bool _ithit = false;
        _ithit = iters > maxit;
        Rprintf("Max iters hit %s\n", _ithit ? "true" : "false");

        meet_cond = ((_thconv && _siconv) || _ithit);

        Rprintf("All convergence %s\n", meet_cond ? "true" : "false");
        converged = _thconv && _siconv;
    }

    // It will be expensive to sequentially grow a list here - do we even need to return
    // all of the intermediate results??
    outlist = List::create(_["FE"]=curr_beta, _["RE"]=curr_u, _["Sigma"]=curr_sigma,
                           _["converged"]=converged, _["Iters"]=iters, _["Dispersion"]=curr_disp,
                           _["Hessian"]=information_sigma);

    return outlist;
}


