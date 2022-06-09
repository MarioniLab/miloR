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
//' @param Z mat - n X (q + g) sparse matrix that maps random effect variable levels to
//' observations - augmented by the n X g genotype matrix
//' @param X mat - n X m sparse matrix that maps fixed effect variables to
//' observations
//' @param K mat - n X n matrix containing genetic relationships between observations
//' @param muvec vec vector of estimated phenotype means
//' @param curr_theta vec vector of initial parameter estimates
//' @param curr_beta vec vector of initial beta estimates
//' @param curr_u vec of initial u estimates
//' @param curr_sigma vec of initial sigma estimates
//' @param curr_G mat c X c matrix of variance components
//' @param y vec of observed counts
//' @param u_indices List a List, each element contains the indices of Z relevant
//' to each RE and all its levels
//' @param theta_conv double Convergence tolerance for paramter estimates
//' @param rlevels List containing mapping of RE variables to individual
//' levels
//' @param curr_disp double Dispersion parameter estimate
//' @param REML bool - use REML for variance component estimation
//' @param maxit int maximum number of iterations if theta_conv is FALSE
//' @param offsets vector of offsets to include in the linear predictor
// [[Rcpp::export]]
List fitGeneticPLGlmm(const arma::mat& Z, const arma::mat& X, const arma::mat& K,
                      arma::vec muvec, arma::vec offsets, arma::vec curr_beta,
                      arma::vec curr_theta, arma::vec curr_u, arma::vec curr_sigma,
                      arma::mat curr_G, const arma::vec& y, List u_indices,
                      double theta_conv,
                      const List& rlevels, double curr_disp, const bool& REML, const int& maxit){

    // declare all variables
    List outlist(11);
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
    V_partial = pseudovarPartial_G(Z, K, u_indices);
    // compute outside the loop
    List VP_partial(c);

    arma::vec score_sigma(c);
    arma::mat information_sigma(c, c);
    arma::vec sigma_update(c);
    arma::vec sigma_diff(sigma_update.size());
    sigma_diff.zeros();

    arma::mat G_inv(stot, stot);

    arma::vec theta_update(m+stot);
    arma::vec theta_diff(theta_update.size());
    theta_diff.zeros();

    List conv_list(maxit+1);

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

    // we only need to invert the Kinship once
    unsigned long _kn = K.n_cols;
    arma::mat Kinv(_kn, _kn);
    // check this isn't singular first (why would it be??)
    double _rcond = arma::rcond(K);
    bool is_singular;
    is_singular = _rcond < 1e-9;

    // check for singular condition
    if(is_singular){
        Rcpp::stop("Kinship matrix is singular");
    }

    Kinv = arma::inv(K); // this could be very slow

    bool converged = false;
    while(!meet_cond){
        D.diag() = muvec; // data space
        Dinv = D.i();
        y_star = computeYStar(X, curr_beta, Z, Dinv, curr_u, y); // data space
        Vmu = computeVmu(muvec, curr_disp);
        W = computeW(curr_disp, Dinv, Vmu);
        Winv = W.i();
        V_star = computeVStar(Z, curr_G, W); // K is implicitly included in curr_G
        V_star_inv = invertPseudoVar(Winv, curr_G, Z);

        if(REML){
            P = computePREML(V_star_inv, X);
            // take the partial derivative outside the while loop, just keep the P*\dVar\dSigma
            VP_partial = pseudovarPartial_P(V_partial, P);

            score_sigma = sigmaScoreREML_arma(VP_partial, y_star, P);
            information_sigma = sigmaInfoREML_arma(VP_partial, P);
        } else{
            List VP_partial = V_partial;
            score_sigma = sigmaScore(y_star, curr_beta, X, VP_partial, V_star_inv);
            information_sigma = sigmaInformation(V_star_inv, VP_partial);
        };

        // information_sigma.brief_print("Fisher\n");
        // double _fishcond = arma::rcond(information_sigma);
        // Rcout << _fishcond << std::endl;

        sigma_update = fisherScore(information_sigma, score_sigma, curr_sigma);
        sigma_diff = abs(sigma_update - curr_sigma); // needs to be an unsigned real value

        // update sigma, G, and G_inv
        curr_sigma = sigma_update;
        curr_G = initialiseG_G(u_indices, curr_sigma, K);
        G_inv = invGmat_G(u_indices, curr_sigma, Kinv);

        // Next, solve pseudo-likelihood GLMM equations to compute solutions for beta and u
        // compute the coefficient matrix
        coeff_mat = coeffMatrix(X, Winv, Z, G_inv); //model space
        // coeff_mat.brief_print("Hessian\n");
        // double _rcond = arma::rcond(coeff_mat);
        // Rcout << _rcond << std::endl;

        theta_update = solveEquations(stot, m, Winv, Z.t(), X.t(), coeff_mat, curr_beta, curr_u, y_star); //model space
        theta_diff = abs(theta_update - curr_theta);

        curr_theta = theta_update; //model space
        curr_beta = curr_theta.elem(beta_ix); //model space
        curr_u = curr_theta.elem(u_ix); //model space
        // curr_u.brief_print("u\n");

        muvec = exp(offsets + (X * curr_beta) + (Z * curr_u)); // data space
        iters++;

        bool _thconv = false;
        _thconv = all(theta_diff < theta_conv);

        bool _siconv = false;
        _siconv = all(sigma_diff < theta_conv);

        bool _ithit = false;
        _ithit = iters > maxit;

        meet_cond = ((_thconv && _siconv) || _ithit);
        converged = _thconv && _siconv;
        List this_conv(5);
        this_conv = List::create(_["ThetaDiff"]=theta_diff, _["SigmaDiff"]=sigma_diff, _["beta"]=curr_beta,
                                 _["u"]=curr_u, _["sigma"]=curr_sigma);
        conv_list(iters-1) = this_conv;
    }

    // inference
    arma::vec se(computeSE(m, stot, coeff_mat));
    arma::vec tscores(computeTScore(curr_beta, se));
    arma::mat vcov(varCovar(VP_partial, c)); // DF calculation is done in R, but needs this

    outlist = List::create(_["FE"]=curr_beta, _["RE"]=curr_u, _["Sigma"]=curr_sigma,
                           _["converged"]=converged, _["Iters"]=iters, _["Dispersion"]=curr_disp,
                           _["Hessian"]=information_sigma, _["SE"]=se, _["t"]=tscores,
                           _["COEFF"]=coeff_mat, _["P"]=P, _["Vpartial"]=VP_partial, _["Ginv"]=G_inv,
                           _["Vsinv"]=V_star_inv, _["Winv"]=Winv, _["VCOV"]=vcov, _["CONVLIST"]=conv_list);

    return outlist;
}


