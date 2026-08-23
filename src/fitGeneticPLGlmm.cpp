#include<RcppArmadillo.h>
#include<cmath>
#include<string>
// [[Rcpp::depends(RcppArmadillo)]]
#include "paramEst.h"
#include "computeMatrices.h"
#include "invertPseudoVar.h"
#include "pseudovarPartial.h"
#include "multiP.h"
#include "inference.h"
#include "utils.h"
using namespace Rcpp;


//' GLMM parameter estimation using pseudo-likelihood with a custom covariance matrix
//'
//' Iteratively estimate GLMM fixed and random effect parameters, and variance
//' component parameters using Fisher scoring based on the Pseudo-likelihood
//' approximation to a Normal loglihood. This function incorporates a user-defined
//' covariance matrix, e.g. a kinship matrix for genetic analyses.
//'
//' @param Z mat - sparse matrix that maps random effect variable levels to
//' observations
//' @param X mat - sparse matrix that maps fixed effect variables to
//' observations
//' @param K mat - sparse matrix that defines the known covariance patterns between
//' individual observations. For example, a kinship matrix will then adjust for the
//' known/estimated genetic relationships between observations.
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
//' @param vardist string which variance form to use NB = negative binomial, P=Poisson [not yet implemented]/
//'
//' @details Fit a NB-GLMM to the counts provided in \emph{y}. The model uses an iterative approach that
//' switches between the joint fixed and random effect parameter inference, and the variance component
//' estimation. A pseudo-likelihood approach is adopted to minimise the log-likelihood of the model
//' given the parameter estimates. The fixed and random effect parameters are estimated using
//' Hendersons mixed model equations, and the variance component parameters are then estimated with
//' the specified solver, i.e. Fisher scoring, Haseman-Elston or constrained Haseman-Elston regression. As
//' the domain of the variance components is [0, +\code{Inf}], any negative variance component estimates will
//' trigger the switch to the HE-NNLS solver until the model converges.
//'
//' @return A \code{list} containing the following elements (note: return types are dictated by Rcpp, so the R
//' types are described here):
//' \describe{
//' \item{\code{FE}:}{\code{numeric} vector of fixed effect parameter estimates.}
//' \item{\code{RE}:}{\code{list} of the same length as the number of random effect variables. Each slot contains the best
//' linear unbiased predictors (BLUPs) for the levels of the corresponding RE variable.}
//' \item{\code{Sigma:}}{\code{numeric} vector of variance component estimates, 1 per random effect variable. For this model the
//' last variance component corresponds to the input \emph{K} matrix.}
//' \item{\code{converged:}}{\code{logical} scalar of whether the model has reached the convergence tolerance or not.}
//' \item{\code{Iters:}}{\code{numeric} scalar with the number of iterations that the model ran for. Is strictly <= \code{max.iter}.}
//' \item{\code{Dispersion:}}{\code{numeric} scalar of the dispersion estimate computed off-line}
//' \item{\code{Hessian:}}{\code{matrix} of 2nd derivative elements from the fixed and random effect parameter inference.}
//' \item{\code{SE:}}{\code{matrix} of standard error estimates, derived from the hessian, i.e. the square roots of the diagonal elements.}
//' \item{\code{t:}}{\code{numeric} vector containing the compute t-score for each fixed effect variable.}
//' \item{\code{COEFF:}}{\code{matrix} containing the coefficient matrix from the mixed model equations.}
//' \item{\code{P:}}{\code{matrix} containing the elements of the REML projection matrix.}
//' \item{\code{Vpartial:}}{\code{list} containing the partial derivatives of the (pseudo)variance matrix with respect to each variance
//' component.}
//' \item{\code{Ginv:}}{\code{matrix} of the inverse variance components broadcast to the full Z matrix.}
//' \item{\code{Vsinv:}}{\code{matrix} of the inverse pseudovariance.}
//' \item{\code{Winv:}}{\code{matrix} of the inverse elements of W = D^-1 V D^-1}
//' \item{\code{VCOV:}}{\code{matrix} of the variance-covariance for all model fixed and random effect variable parameter estimates.
//' This is required to compute the degrees of freedom for the fixed effect parameter inference.}
//' \item{\code{CONVLIST:}}{\code{list} of \code{list} containing the parameter estimates and differences between current and previous
//' iteration estimates at each model iteration. These are included for each fixed effect, random effect and variance component parameter.
//' The list elements for each iteration are: \emph{ThetaDiff}, \emph{SigmaDiff}, \emph{beta}, \emph{u}, \emph{sigma}.}
//' }
//'
//' @author Mike Morgan
//'
//' @examples
//' NULL
//'
//' @name fitGeneticPLGlmm
//'
// [[Rcpp::export]]
List fitGeneticPLGlmm(const arma::mat& Z, const arma::mat& X, const arma::mat& K,
                      arma::vec muvec, arma::vec offsets, arma::vec curr_beta,
                      arma::vec curr_theta, arma::vec curr_u, arma::vec curr_sigma,
                      arma::mat curr_G, const arma::vec& y, List u_indices,
                      double theta_conv,
                      const List& rlevels, double curr_disp, const bool& REML, const int& maxit,
                      std::string solver,
                      std::string vardist){

    // no guarantee that Pi exists before C++ 20(?!?!?!)
    constexpr double pi = 3.14159265358979323846;

    // declare all variables
    List outlist(14);
    int iters=0;
    int stot = Z.n_cols;
    const int& c = curr_sigma.size();
    const int& m = X.n_cols;
    const int& n = X.n_rows;
    bool meet_cond = false;
    double constval = 1e-8; // value at which to constrain values
    double _intercept = constval; // intercept for HE regression
    double max_disp = 1e4; // ceiling on the dispersion search, as in fitGeneticNullGlmm
    double delta_up = std::min(max_disp, 2.0 * curr_disp);
    double delta_lo = 1e-2;
    double update_disp = 0.0;
    double disp_diff = 0.0;

    // setup matrices
    arma::mat D(n, n);
    D.zeros();
    arma::mat Dinv(n, n);
    Dinv.zeros();

    arma::vec y_star(n);
    arma::vec y_star_c(n); // working response with the offset removed

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
    List precomp_list(2);
    List pzzp_list(c); // P * Z(j) * Z(j)^T * P^T

    arma::vec score_sigma(c);
    arma::mat information_sigma(c, c);
    arma::vec sigma_update(c);
    arma::vec _sigma_update(c+1);
    arma::vec sigma_diff(sigma_update.size());
    sigma_diff.zeros();

    arma::mat G_inv(stot, stot);
    arma::mat Zstar(n, stot);
    arma::mat Gfill(stot, stot);

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

    // The genetic random effect is carried in the whitened parameterisation:
    // Z_g is the Cholesky factor L of the relatedness matrix, set in fitGLMM, so
    // Z_g Z_g' = L L' = K and the corresponding block of G is simply sigma_g I.
    // G used to hold sigma_g K in that block instead, which - against a Z_g that
    // is already L - built sigma_g L K L' into the pseudo-variance and applied K
    // twice. Every partial derivative of V* with respect to sigma_g in this file
    // is K, so the score and the information described a model the pseudo-
    // variance did not encode. In the whitened form G and G_inv are diagonal, so
    // no inverse of K is needed here at all.
    if(K.n_cols != n || K.n_rows != n){
        stop("Covariance matrix and design matrix dimensions do not match");
    }

    arma::uvec _gidx = u_indices[c - 1];
    if(_gidx.n_elem != static_cast<unsigned int>(n)){
        stop("RE indices and dimensions of covariance do not match");
    }

    bool converged = false;
    bool _phi_est = true; // control if we re-estimate phi or not

    // initial optimisation of dispersion
    // switch this to a golden-section search?
    update_disp = phiGoldenSearch(curr_disp, delta_lo, delta_up, c,
                                  muvec, G_inv, pi,
                                  curr_u, curr_sigma, y);
    disp_diff = abs(curr_disp - update_disp);
    // curr_disp = update_disp;
    // make the upper and lower bounds based on the current value,
    // but 0 < lo < up < 1.0
    // Re-bracket around the value the search just returned, on both sides of
    // it. The bracket used to be [phi/2, phi], whose upper end is the current
    // estimate itself, so the golden section search could only ever return a
    // smaller value: the dispersion halved every iteration regardless of data.
    delta_lo = std::max(1e-2, update_disp * 0.5);
    delta_up = std::min(max_disp, std::max(update_disp * 2.0, delta_lo + 1e-2));

    while(!meet_cond){
        curr_disp = update_disp;
        D.diag() = muvec; // data space
        Dinv = D.i();
        y_star = computeYStar(X, curr_beta, Z, Dinv, curr_u, y, offsets); // data space
        // The offset belongs to the linear predictor but is not a column of X,
        // so it has to come off the working response before it is used: P only
        // annihilates the column space of X, and the mixed model equations
        // otherwise let the intercept absorb the offset and feed it back into
        // eta a second time on the next iteration.
        y_star_c = y_star - offsets;

        Vmu = computeVmu(muvec, curr_disp, vardist);
        W = computeW(curr_disp, Dinv, vardist);
        Winv = W.i();

        // pre-compute matrics: X^T * W^-1, Z^T * W^-1
        arma::mat xTwinv = X.t() * Winv;
        arma::mat zTwin = Z.t() * Winv;

        V_star = computeVStar(Z, curr_G, W); // K is implicitly included in curr_G
        V_star_inv = invertPseudoVar(Winv, curr_G, Z, zTwin);

        if(REML){
            P = computePREML(V_star_inv, X);
            // take the partial derivative outside the while loop, just keep the P*\dVar\dSigma
        } else{
            P = arma::eye(n, n);
        };

        // pre-compute matrics: P*Z
        arma::mat PZ = P * Z;

        // pre-compute P*Z(j) * Z(j)^T
        precomp_list = computePZList_G(u_indices, PZ, P, Z, solver, K);

        // choose between HE regression and Fisher scoring for variance components
        // sigma_update is always 1 element longer than the others with HE, but we need to keep track of this
        if(solver == "HE"){
            // try Haseman-Elston regression instead of Fisher scoring
            sigma_update = estHasemanElstonGenetic(Z, P, PZ, u_indices, y_star_c, K, W);
        } else if (solver == "HE-NNLS"){
            // for the first iteration use the current non-zero estimate
            arma::dvec _curr_sigma(c+1, arma::fill::zeros);

            if(REML){
                _sigma_update = estHasemanElstonConstrainedGenetic(Z, P, PZ, u_indices, y_star_c, K, _curr_sigma, iters, W);
                _intercept = _sigma_update[0];
                sigma_update = _sigma_update.tail(c);
            } else{
                _sigma_update = estHasemanElstonConstrainedGeneticML(Z, u_indices, y_star_c, K, _curr_sigma, iters, W);
                _intercept = _sigma_update[0];
                sigma_update = _sigma_update.tail(c);
            }

            // set 0 values to minval to prevent 0 denominators later
            if(any(sigma_update == 0.0)){
                for(int i=0; i<c; i++){
                    if(sigma_update[i] <= 0.0){
                        sigma_update[i] = constval;
                    }
                }
            }

        }else if(solver == "Fisher"){
            if(REML){
                // Score and information both in the P basis - the score used to
                // mix P for its trace term with Vstar^-1 for its quadratic form,
                // which is not the derivative of any objective.
                VP_partial = precomp_list["PZZt"];

                score_sigma = sigmaScoreREML_arma(VP_partial, y_star_c, P,
                                                  curr_beta, X);
                information_sigma = sigmaInfoREML_arma(VP_partial, P);
            } else{
                // this used to redeclare VP_partial, shadowing the outer list, so
                // the Vpartial returned to R was empty for ML with Fisher scoring
                // and varCovar() ran on nothing
                VP_partial = V_partial;
                score_sigma = sigmaScore(y_star_c, curr_beta, X, VP_partial, V_star_inv);
                information_sigma = sigmaInformation(V_star_inv, VP_partial);
            }
            sigma_update = fisherScore(information_sigma, score_sigma, curr_sigma);

            // The domain of the variance components is [0, Inf). Retreat along
            // the ascent direction by step-halving rather than abandoning Fisher
            // scoring the moment a full step oversteps the boundary - this keeps
            // the search direction and lets a component recover on a later
            // iteration if the data support it.
            arma::vec fisher_step = sigma_update - curr_sigma;
            int halvings = 0;
            while(arma::any((curr_sigma + fisher_step) <= 0.0) && halvings < 30){
                fisher_step *= 0.5;
                halvings++;
            }
            sigma_update = curr_sigma + fisher_step;

            // final guard for non-finite or still non-positive components
            for(int i=0; i < c; i++){
                if(!std::isfinite(sigma_update[i]) || sigma_update[i] <= 0.0){
                    sigma_update[i] = constval;
                }
            }
        }

        // if we have negative sigmas then we need to switch solver
        if(any(sigma_update < 0.0)){
            warning("Negative variance components - re-running with NNLS");
            solver = "HE-NNLS";
            // for the first iteration use the current non-zero estimate
            arma::dvec _curr_sigma(c+1, arma::fill::zeros);

            if(REML){
                _sigma_update = estHasemanElstonConstrainedGenetic(Z, P, PZ, u_indices, y_star_c, K, _curr_sigma, iters, W);
                _intercept = _sigma_update[0];
                sigma_update = _sigma_update.tail(c);
            } else{
                _sigma_update = estHasemanElstonConstrainedGeneticML(Z, u_indices, y_star_c, K, _curr_sigma, iters, W);
                _intercept = _sigma_update[0];
                sigma_update = _sigma_update.tail(c);
            }

            // set 0 values to minval to prevent 0 denominators later
            if(any(sigma_update == 0.0)){
                for(int i=0; i<c; i++){
                    if(sigma_update[i] <= 0.0){
                        sigma_update[i] = constval;
                    }
                }
            }
        }

        sigma_diff = abs(sigma_update - curr_sigma); // needs to be an unsigned real value

        // update sigma, G, and G_inv
        curr_sigma = sigma_update;
        curr_G = initialiseG(u_indices, curr_sigma);
        G_inv = invGmat(u_indices, curr_sigma);

        // Update the dispersion with the new variances.
        //
        // This used to be skipped once two consecutive searches agreed to
        // within 1e-2, but that is a convergence test on the wrong quantity:
        // the dispersion is estimated from the NB likelihood at the current mu,
        // and mu is still moving in the early iterations. Two early searches
        // agreeing switched the estimator off permanently and left the
        // dispersion pinned to a mu the model had already left behind, which
        // then propagates into W and starves the variance components estimated
        // on the same pseudo-variance.
        update_disp = phiGoldenSearch(curr_disp, delta_lo, delta_up, c,
                                      muvec, G_inv, pi,
                                      curr_u, curr_sigma, y);

        // bracket the current estimate on both sides - see above
        delta_lo = std::max(1e-2, update_disp * 0.5);
        delta_up = std::min(max_disp, std::max(update_disp * 2.0, delta_lo + 1e-2));

        // update_disp = phiMME(y_star, curr_sigma);
        disp_diff = abs(curr_disp - update_disp);

        // Next, solve pseudo-likelihood GLMM equations to compute solutions for beta and u
        // compute the coefficient matrix
        coeff_mat = coeffMatrix(X, xTwinv, zTwin, Z, G_inv); //model space
        theta_update = solveEquations(stot, m, zTwin, xTwinv, coeff_mat, curr_beta, curr_u, y_star_c); //model space

        LogicalVector _check_theta = check_na_arma_numeric(theta_update);
        bool _any_ystar_na = any(_check_theta).is_true(); // .is_true required for proper type casting to bool

        if(_any_ystar_na){
            List this_conv(5);
            this_conv = List::create(_["ThetaDiff"]=theta_diff, _["SigmaDiff"]=sigma_diff, _["beta"]=curr_beta,
                                     _["u"]=curr_u, _["sigma"]=curr_sigma);
            conv_list(iters-1) = this_conv;
            warning("NaN in theta update");
            break;
        }

        theta_diff = abs(theta_update - curr_theta);

        curr_theta = theta_update; //model space
        curr_beta = curr_theta.elem(beta_ix); //model space
        curr_u = curr_theta.elem(u_ix); //model space

        muvec = exp(offsets + (X * curr_beta) + (Z * curr_u)); // data space
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

        // compute final loglihood
        // make non-broadcast G matrix
        arma::mat littleG(c, c, arma::fill::zeros);

        for(int i=0; i<c; i++){
            littleG(i, i) = curr_sigma(i);
        }
        double loglihood = nbLogLik(muvec, curr_disp, y) - normLogLik(c, G_inv, littleG, curr_u, pi);

        List this_conv(8);
        this_conv = List::create(_["ThetaDiff"]=theta_diff, _["SigmaDiff"]=sigma_diff, _["beta"]=curr_beta,
                                 _["u"]=curr_u, _["sigma"]=curr_sigma, _["disp"]=curr_disp, _["PhiDiff"]=disp_diff,
                                 _["LOGLIHOOD"]=loglihood);
        conv_list(iters-1) = this_conv;
    }

    // inference
    arma::vec se(computeSE(m, stot, coeff_mat));
    arma::vec tscores(computeTScore(curr_beta, se));

    // VP_partial is empty for HE/HE-NNLS so needs to be populated
    if(solver != "Fisher"){
        VP_partial = precomp_list["PZZt"];
    }

    arma::mat vcov(varCovar(VP_partial, c)); // DF calculation is done in R, but needs this

    // compute the variance of the pseudovariable
    double pseduo_var = arma::var(y_star);

    // compute final loglihood
    // make non-broadcast G matrix
    arma::mat littleG(c, c, arma::fill::zeros);

    for(int i=0; i<c; i++){
        littleG(i, i) = curr_sigma(i);
    }
    double loglihood = nbLogLik(muvec, curr_disp, y) - normLogLik(c, G_inv, littleG, curr_u, pi);

    outlist = List::create(_["FE"]=curr_beta, _["RE"]=curr_u, _["Sigma"]=curr_sigma,
                           _["converged"]=converged, _["Iters"]=iters, _["Dispersion"]=curr_disp,
                           _["Hessian"]=information_sigma, _["SE"]=se, _["t"]=tscores, _["PSVAR"]=pseduo_var,
                           _["COEFF"]=coeff_mat, _["P"]=P, _["Vpartial"]=VP_partial, _["Ginv"]=G_inv,
                           _["Vsinv"]=V_star_inv, _["Winv"]=Winv, _["VCOV"]=vcov, _["LOGLIHOOD"]=loglihood,
                            _["CONVLIST"]=conv_list);

    return outlist;
}


