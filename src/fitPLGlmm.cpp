#include<RcppArmadillo.h>
#include<Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "paramEst.h"
#include "computeMatrices.h"
#include "invertPseudoVar.h"
#include "pseudovarPartial.h"
#include "multiP.h"
#include "inference.h"
#include "utils.h"
using namespace Rcpp;

//' GLMM parameter estimation using pseudo-likelihood
//'
//' Iteratively estimate GLMM fixed and random effect parameters, and variance
//' component parameters using Fisher scoring based on the Pseudo-likelihood
//' approximation to a Normal loglihood.
//' @param Z mat - sparse matrix that maps random effect variable levels to
//' observations
//' @param X mat - sparse matrix that maps fixed effect variables to
//' observations
//' @param muvec vec vector of estimated phenotype means
//' @param offsets vec vector of model offsets
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
//' @param solver string which solver to use - either HE (Haseman-Elston regression) or Fisher scoring
// [[Rcpp::export]]
List fitPLGlmm(const arma::mat& Z, const arma::mat& X, arma::vec muvec,
               arma::vec offsets, arma::vec curr_beta,
               arma::vec curr_theta, arma::vec curr_u, arma::vec curr_sigma,
               arma::mat curr_G, const arma::vec& y, List u_indices,
               double theta_conv,
               const List& rlevels, double curr_disp, const bool& REML, const int& maxit,
               std::string solver){

    // declare all variables
    List outlist(10);
    int iters=0;
    int stot = Z.n_cols;
    const int& c = curr_sigma.size();
    const int& m = X.n_cols;
    const int& n = X.n_rows;
    bool meet_cond = false;
    double constval = 1e-8; // value at which to constrain values
    double _intercept = constval; // intercept for HE regression

    // setup matrices
    arma::mat D(n, n);
    D.zeros();
    arma::mat Dinv(n, n);
    Dinv.zeros();

    arma::vec y_star(n);

    arma::mat Vmu(n, n);
    Vmu.zeros();
    arma::mat W(n, n);
    W.zeros();
    arma::mat Winv(n, n);
    Winv.zeros();

    arma::mat V_star(n, n);
    V_star.zeros();
    arma::mat V_star_inv(n, n);
    V_star_inv.zeros();
    arma::mat P(n, n);
    P.zeros();

    arma::mat coeff_mat(m+c, m+c);
    coeff_mat.zeros();
    List V_partial(c);
    V_partial = pseudovarPartial_C(Z, u_indices);
    // compute outside the loop
    List VP_partial(c);

    arma::vec score_sigma(c);
    arma::mat information_sigma(c, c);
    information_sigma.zeros();
    arma::vec sigma_update(c);
    arma::vec sigma_diff(sigma_update.size());
    sigma_diff.zeros();

    arma::mat G_inv(stot, stot);
    G_inv.zeros();

    arma::vec theta_update(m+stot);
    arma::vec theta_diff(theta_update.size());
    theta_diff.zeros();

    List conv_list(maxit+1);

    // create a uvec of sigma indices
    arma::uvec _sigma_index(c);
    for(int i=0; i < c; i++){
        _sigma_index[i] = i+1;
    }

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

        // check for all zero eigen values
        arma::cx_vec d_eigenval = arma::eig_gen(D); // this needs to handle complex values
        LogicalVector _check_zero = check_zero_arma_complex(d_eigenval);
        bool _all_zero = any(_check_zero).is_true();

        if(_all_zero){
            stop("Zero eigenvalues in D - do you have collinear variables?");
        }

        Dinv = D.i();
        y_star = computeYStar(X, curr_beta, Z, Dinv, curr_u, y, offsets);
        Vmu = computeVmu(muvec, curr_disp);
        W = computeW(curr_disp, Dinv, Vmu);
        Winv = W.i();
        V_star = computeVStar(Z, curr_G, W);
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

        // choose between HE regression and Fisher scoring for variance components
        // would a hybrid approach work here? If any HE estimates are zero switch
        // to NNLS using these as the initial estimates?
        if(solver == "HE"){
            // try Haseman-Elston regression instead of Fisher scoring
            sigma_update = estHasemanElston(Z, P, u_indices, y_star);
        } else if(solver == "HE-NNLS"){
            // for the first iteration use the current non-zero estimate
            arma::dvec _curr_sigma(c+1, arma::fill::zeros);
            _curr_sigma.fill(constval);

            // if these are all zero then it can only be that they are initial estimates
            if(iters > 0){
                _curr_sigma[0] = _intercept;
                _curr_sigma.elem(_sigma_index) = curr_sigma; // is this valid to set elements like this?
            }
            sigma_update = estHasemanElstonConstrained(Z, P, u_indices, y_star, _curr_sigma, iters);
        }else if(solver == "Fisher"){
            sigma_update = fisherScore(information_sigma, score_sigma, curr_sigma);
        }

        // if we have negative sigmas then we need to switch solver
        if(any(sigma_update < 0.0)){
            warning("Negative variance components - re-running with NNLS");
            solver = "HE-NNLS";
            // for the first iteration use the current non-zero estimate
            arma::dvec _curr_sigma(c+1, arma::fill::zeros);
            _curr_sigma.fill(constval);

            // if these are all zero then it can only be that they are initial estimates
            if(iters > 0){
                _curr_sigma[0] = _intercept;
                _curr_sigma.elem(_sigma_index) = curr_sigma; // is this valid to set elements like this?
            }
            sigma_update = estHasemanElstonConstrained(Z, P, u_indices, y_star, _curr_sigma, iters);
        }

        // update sigma, G, and G_inv
        curr_sigma = sigma_update;
        curr_G = initialiseG(u_indices, curr_sigma);
        G_inv = invGmat(u_indices, curr_sigma);

        // Next, solve pseudo-likelihood GLMM equations to compute solutions for B and u
        // compute the coefficient matrix
        coeff_mat = coeffMatrix(X, Winv, Z, G_inv);
        theta_update = solveEquations(stot, m, Winv, Z.t(), X.t(), coeff_mat, curr_beta, curr_u, y_star);
        theta_diff = abs(theta_update - curr_theta);

        // inference
        curr_theta = theta_update;
        curr_beta = curr_theta.elem(beta_ix);
        curr_u = curr_theta.elem(u_ix);

        // need to check for infinite and NA values here...
        // muvec = exp(offsets + (X * curr_beta) + (Z * curr_u));
        muvec = exp(offsets + (X * curr_beta) + (Z * curr_u));
        LogicalVector _check_mu = check_na_arma_numeric(muvec);
        bool _any_na = any(_check_mu).is_true(); // .is_true required for proper type casting to bool

        LogicalVector _check_inf = check_inf_arma_numeric(muvec);
        bool _any_inf = any(_check_inf).is_true();

        if(_any_na){
            stop("NA estimates in linear predictor - consider an alternative model");
        }

        if(_any_inf){
            stop("Infinite parameter estimates - consider an alternative model");
        }

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
