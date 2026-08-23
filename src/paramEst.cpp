#define ARMA_WARN_LEVEL 1
#include "paramEst.h"
#include "computeMatrices.h"
#include "utils.h"
#include<RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include<cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// using namespace Rcpp;

// All functions used in parameter estimation

arma::vec sigmaScoreREML_arma (const Rcpp::List& PdV, const arma::vec& ystar,
                               const arma::mat& P, const arma::vec& curr_beta,
                               const arma::mat& X){
    // REML score for the variance components:
    //
    //   dl/dsigma_j = -0.5 tr(P dV/dsigma_j) + 0.5 y*' P dV/dsigma_j P y*
    //
    // Both terms are in the projection basis P. Taking the trace in P and the
    // quadratic form in Vstar^-1 does not differentiate any objective: Vstar^-1
    // - P is positive semi-definite, so the quadratic form is inflated relative
    // to the trace it is meant to balance and the root of the score moves away
    // from the REML solution. Only PdV is needed here, so there is no second
    // basis to pass in and get wrong.
    //
    // ystar must have the offset removed: the offset is part of the linear
    // predictor but not a column of X, so P does not annihilate it.
    const int& c = PdV.size();
    const int& n = X.n_rows;
    arma::vec reml_score(c);
    arma::vec ystarminx(n);
    ystarminx = ystar - (X * curr_beta);

    // P X = 0, so this is P y* regardless of the current beta
    arma::vec Presid = P * ystarminx;

    for(int i=0; i < c; i++){
        const arma::mat& PdVi = PdV(i); // this is P * partial derivative

        double lhs = -0.5 * arma::trace(PdVi);
        // ystarminx' (P dV_i) (P ystarminx) = ystarminx' P dV_i P ystarminx
        double rhs = 0.5 * arma::dot(ystarminx, PdVi * Presid);

        reml_score[i] = lhs + rhs;
    }

    return reml_score;
}


arma::mat sigmaInfoREML_arma (const Rcpp::List& pvstari, const arma::mat& P){
    // REML Fisher/expected information matrix
    // sparsifying this is baaaad for performance
    const int& c = pvstari.size();
    arma::mat sinfo(c, c);

    // this is a symmetric matrix so only need to fill the upper or
    // lower triangle to make it O(n^2/2) rather than O(n^2)
    for(int i=0; i < c; i++){
        const arma::mat& _ipP = pvstari[i]; // this is P * \d Var/ \dsigma

        for(int j=i; j < c; j++){
            const arma::mat& P_jp = pvstari[j]; // this is P * \d Var/ \dsigma

            // Compute trace directly
            double trace = 0.0;
            for(arma::uword k = 0; k < _ipP.n_rows; ++k) {
                trace += arma::dot(_ipP.row(k), P_jp.col(k));
            }

            sinfo(i, j) = 0.5 * trace;
            if(i != j){
                sinfo(j, i) = 0.5 * trace;
            }

        }
    }

    return sinfo;
}


arma::vec sigmaScore (arma::vec ystar, arma::vec beta, arma::mat X, Rcpp::List V_partial, arma::mat V_star_inv){

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
        score[i] = lhs + 0.5 * rhs_mat(0, 0);
    }

    return score;
}


arma::mat sigmaInformation (arma::mat V_star_inv, Rcpp::List V_partial){
    int c = V_partial.size();
    int n = V_star_inv.n_cols;
    arma::mat sinfo = arma::zeros(c, c);

    for(int i=0; i < c; i++){
        arma::mat _ip = V_partial(i);
        for(int j=0; j < c; j++){
            arma::mat _jp = V_partial(j);

            arma::mat _inmat(n, n);
            double _tr = 0.0;
            _inmat = V_star_inv * _ip * V_star_inv * _jp;
            _tr = 0.5 * arma::trace(_inmat);

            sinfo(i, j) = _tr;
        }
    }

    return sinfo;
}


arma::vec fisherScore (const arma::mat& hess, const arma::vec& score_vec, const arma::vec& theta_hat){
    // sequentially update the parameter using the Newton-Raphson algorithm
    // theta ~= theta_hat + hess^-1 * score
    // this needs to be in a direction of descent towards a minimum
    int m = theta_hat.size();
    arma::vec theta(m, arma::fill::zeros);
    arma::mat hessinv(hess.n_rows, hess.n_cols);

    // will need a check here for singular hessians...
    double _rcond = arma::rcond(hess);
    if(_rcond < 1e-9){
        Rcpp::warning("Variance Component Hessian is computationally singular");
        hessinv = arma::pinv(hess);
    } else{
        hessinv = arma::inv(hess); // always use pinv? solve() and inv() are most sensitive than R versions
    }
    theta = theta_hat + (hessinv * score_vec);

    return theta;
}


arma::mat coeffMatrix(const arma::mat& X, const arma::mat& XtWinv, const arma::mat& ZtWinv,
                      const arma::mat& Z, const arma::mat& Ginv){
    // compute the components of the coefficient matrix for the MMEs
    // sparsification _does_ help here, despite the added overhead
    // unsparsify
    int c = Z.n_cols;
    int m = X.n_cols;
    arma::mat lhs(m+c, m+c);

    lhs(arma::span(0, m-1), arma::span(0, m-1)) = XtWinv * X;
    lhs(arma::span(0, m-1), arma::span(m, m+c-1)) = XtWinv * Z;
    lhs(arma::span(m, m+c-1), arma::span(0, m-1)) = ZtWinv * X;
    lhs(arma::span(m, m+c-1), arma::span(m, m+c-1)) = ZtWinv * Z + Ginv;

    return lhs;
}


arma::vec solveEquations (const int& c, const int& m, const arma::mat& ZtWinv, const arma::mat& XtWinv,
                          const arma::mat& coeffmat, const arma::vec& beta, const arma::vec& u, const arma::vec& ystar){
    // solve the mixed model equations
    arma::vec rhs_beta(m);
    arma::vec rhs_u(c);
    arma::mat rhs(m+c, 1);

    arma::vec theta_up(m+c, arma::fill::zeros);
    arma::mat I = arma::eye(arma::size(coeffmat));
    arma::mat lm_eye = arma::eye(arma::size(coeffmat));
    arma::mat _coeff = coeffmat;
    double lambda = 1e-1;
    double lambda_step = 10;
    double _illcond_eps = 1e-6;
    double _lcond_target = 1e-5;

    rhs_beta.col(0) = XtWinv * ystar;
    rhs_u.col(0) = ZtWinv * ystar;

    rhs = arma::join_cols(rhs_beta, rhs_u);

    // need a check for singular hessian here
    double _rcond = arma::rcond(coeffmat);

    try{
        if(_rcond < 1e-9){
            Rcpp::warning("Coefficients Hessian is computationally singular - trying pseudoinverse");
            // poorly conditioned system - switch to regularisation to find a solution?
            // add a small diagonal element to coeffmat
            theta_up = arma::solve(_coeff, rhs, arma::solve_opts::no_approx);
        } else{
            // can we just use solve here instead?
            // if the coefficient matrix is singular then do we resort to pinv?
            theta_up = arma::solve(coeffmat, rhs, arma::solve_opts::fast);
        }
    } catch(std::exception const& ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }
    return theta_up;
}


arma::mat computeZstar(const arma::mat& Z, const arma::vec& curr_sigma, const Rcpp::List& u_indices){
    // use the A = LL^T, ZL = Z* trick similar to lme4 - described in more detail
    // in https://onlinelibrary.wiley.com/doi/epdf/10.1046/j.1439-0388.2002.00327.x
    // this assumes that G is PSD - it will fail if it is not

    // compute the Cholesky of G
    int stot = Z.n_cols;
    int n = Z.n_rows;
    arma::mat G(stot, stot);
    G = initialiseG(u_indices, curr_sigma);

    arma::mat cholG(stot, stot);
    try{
        cholG = arma::chol(G, "lower");
    } catch(std::exception &ex){
        Rcpp::stop("G is not positive (semi) definite - Cholesky failed");
        // forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }

    // construct Z*=ZL
    arma::mat Zstar(n, stot);
    Zstar = Z * cholG;

    return Zstar;
}


arma::vec nnlsSolve(const arma::mat& vecZ, const arma::vec& Y, arma::vec nnls_update, const int& Iters){
    // Lawson and Hanson algorithm for constrained NNLS
    // initial params
    // P is the passive set - i.e. the indices of the non-negative estimates
    // R is the active set - the indices of the constrained estimates, i.e. those held at zero
    // d = 0 is a solution vector
    // w is a vector of the Langrangian multipliers

    // need to implement a check for infinite loop - have a max iters argument?
    double constval = 0.0; // value at which to constrain values
    double EPS = 1e-6;
    unsigned int m = vecZ.n_cols;

    // Lawson and Hanson start from x = 0 with an empty passive set and a full
    // active set. Only a strictly positive component is free; a component sitting
    // at the boundary belongs to the active set. These two tests used to be on
    // the wrong side of zero, so a zero starting vector - which is what both the
    // HE-NNLS solver and the negative-variance fallback hand in - put every index
    // in the passive set and left the active set empty. An empty active set is
    // itself a termination condition, so the loop below was skipped entirely and
    // the zero vector was handed straight back as the variance component estimate.
    arma::ivec P(m, arma::fill::ones); // the indices have to be set to negative values to be empty

    for(int i=0; i < m; i++){
        if(nnls_update[i] <= constval){
            P[i] = -1;
        }
    }

    arma::ivec R(m, arma::fill::ones); // these are the indices of the active set

    for(int i=0; i < m; i++){
        if(nnls_update[i] > constval){
            R[i] = -1;
        }
    }

    arma::dvec w(m);
    w = vecZ.t() * (Y - vecZ*nnls_update); //Lagrangian multipliers

    bool check_R;
    bool check_wR;
    bool check_wP;
    bool check_conditions;
    check_R = all(R < 0);
    check_wR = all(w.elem(find(R > 0)) <= constval); // these should be 0 or negative
    check_wP = all(abs(w.elem(find(P > 0)) - constval) < EPS); // these should be ~0
    check_conditions = check_R || (check_wR && check_wP);

    double max_w = 0.0;

    // Lawson and Hanson terminates in a finite number of steps in exact
    // arithmetic, but the loop below has no bound of its own and this routine
    // used to return before entering it, so the absence of one was never
    // exercised. Cap it: each pass frees one index, so a few sweeps of the
    // parameter vector is a generous bound.
    const unsigned int max_outer = 3 * m + 10;
    unsigned int outer = 0;

    while(!check_conditions && outer < max_outer){ // set a tolerance here to check for positive Langrangian
        outer++;
        arma::uvec active_idx = find(R > 0);
        max_w = max(w.elem(active_idx)); // what if several are identical?
        unsigned int max_j = active_idx(0);

        // get the index of the maximal Langrangian multiplier. The search has to
        // be restricted to the active set: taken over every index it can land on
        // a passive one whose multiplier happens to sit within EPS of the
        // maximum, which flips the wrong index in and out of the passive set and
        // leaves the loop cycling.
        for(arma::uword k=0; k < active_idx.n_elem; k++){
            if(abs(w[active_idx(k)] - max_w) <= EPS){
                max_j = active_idx(k);
            }
        }

        P[max_j] = -1 * P[max_j]; // this should reverse the sign
        R[max_j] = -1 * R[max_j];

        // find the elements >= 0
        arma::uvec select_P = find(P > 0);
        arma::vec s_all(m);
        s_all.fill(constval);
        s_all.elem(select_P) = arma::solve(vecZ.cols(select_P), Y, arma::solve_opts::fast); // is this faster?

        double min_sp = s_all.elem(select_P).min();
        double alpha = 1.0; // the step size
        unsigned int inner = 0;

        while(min_sp < 0 && inner < max_outer){
            inner++;
            // recompute selection of P element restricted to negative estimates in S
            arma::uvec select_sP = find(P > 0 && s_all < 0);
            arma::vec diffVec(select_sP.size());

            // zero divisions create problems here
            // need a check for inf, -inf and nan
            diffVec = nnls_update.elem(select_sP)/(nnls_update.elem(select_sP) - s_all.elem(select_sP));
            alpha = diffVec.min();

            // update estimates
            nnls_update = nnls_update + (alpha * (s_all - nnls_update));

            // switch any zeros to P and negatives to R
            arma::uvec _isswitch(m, arma::fill::ones);

            for(int i=0; i < m; i++){
                if(nnls_update[i] <= constval){
                    _isswitch[i] = 0;
                }
            }

            arma::uvec switch_j = find(_isswitch == 0);

            P.elem(switch_j) = -1 * P.elem(switch_j);
            R.elem(switch_j) = -1 * R.elem(switch_j);
            select_P = find(P > 0);

            // update the elements of S
            s_all.fill(constval);
            s_all.elem(select_P) = arma::solve(vecZ.cols(select_P), Y, arma::solve_opts::fast);
            // select_P needs to be non-empty
            // if P is empty then we don't have any non-negative elements.
            // if select_P is empty then end

            if(select_P.n_rows == 0){
                min_sp = constval;
            } else{
                min_sp = s_all.elem(select_P).min();
            }
        }

        nnls_update = s_all;
        w = vecZ.t() * (Y - vecZ * nnls_update);
        check_R = all(R < 0);
        check_wR = all(w.elem(find(R > 0)) <= constval); // these should be 0 or negative
        check_wP = all(abs(abs(w.elem(find(P > 0))) - constval) < EPS); // these should be ~0
        check_conditions = check_R || (check_wR && check_wP);
    }

    // at convergence w^R < 0 and w^P = 0
    return nnls_update;
}


arma::vec fastNnlsSolve(const arma::mat& vecZ, const arma::vec& Y){
    // This uses the Fast Approximate Solution Trajectory NNLS
    // from https://stackoverflow.com/questions/58006606/rcpp-implementation-of-fast-nonnegative-least-squares

    // assumes the problem is not ill-conditioned - how realistic is this?
    // vecZ is almost _never_ square
    int m = vecZ.n_rows;
    arma::mat A(m, m);
    A = vecZ * vecZ.t();
    arma::vec nnls_update = arma::solve(A, Y,
                                        arma::solve_opts::likely_sympd + arma::solve_opts::fast);

    while(any(nnls_update < 0)){
        // define the feasible set P as all estimates > 0
        arma::uvec nonzero = find(nnls_update > 0);

        // reset estimates
        nnls_update.zeros();

        // now solve the OLS problem for params in the feasible set
        nnls_update(nonzero) = arma::solve(vecZ.submat(nonzero, nonzero), Y.elem(nonzero),
                    arma::solve_opts::likely_sympd + arma::solve_opts::fast);
    }
    return nnls_update;
}


double phiLineSearch(double disp, double lower, double upper, const int& c,
                     const arma::vec& mu, const arma::mat& Ginv, double pi,
                     const arma::vec& curr_u, const arma::vec& sigma,
                     const arma::vec& y){
    // perform a bisection search for dispersion
    // evaluate the loglihood at each bound
    arma::mat littleG(c, c, arma::fill::zeros);

    for(int i=0; i<c; i++){
        littleG(i, i) = sigma(i);
    }

    double normlihood = normLogLik(c, Ginv, littleG, curr_u, pi);
    double half_logli = nbLogLik(mu, disp/2.0, y) - normlihood;
    double up_logli = nbLogLik(mu, lower, y) - normlihood;
    double lo_logli = nbLogLik(mu, upper, y) - normlihood;
    double new_disp = 0.0;

    bool comp_vals = false;

    if(lo_logli < up_logli){
        new_disp = lower;
    } else{
        new_disp = upper;
    }

    return new_disp;
}


double phiGoldenSearch(double disp, double lower, double upper, const int& c,
                       const arma::vec& mu, const arma::mat& Ginv, double pi,
                       const arma::vec& curr_u, const arma::vec& sigma,
                       const arma::vec& y){

    // perform a golden section search
    // function inspired by https://drlvk.github.io/nm/section-golden-section.html
    // this only requires a single call to this function

    // Tolerance on the width of the bracket, relative to its scale. This was a
    // flat 1e-2, which leaves the returned value uncertain at the 1e-3 level.
    // The dispersion is re-estimated on every iteration, so that uncertainty
    // reappears as jitter in W and prevents the outer fixed point from ever
    // meeting a parameter tolerance tighter than it. The search is O(n) per
    // evaluation against the O(n^3) work in the same iteration, so resolving it
    // properly costs roughly twenty extra evaluations and nothing else.
    double tol = 1e-8 * std::max(1.0, upper);
    double r = (3 - std::sqrt(5))/2.0;
    double pc = lower + r*(upper - lower);
    double pd = upper - r*(upper - lower);

    // make non-broadcast G matrix
    arma::mat littleG(c, c, arma::fill::zeros);

    for(int i=0; i<c; i++){
        littleG(i, i) = sigma(i);
    }

    // evaluate normal likelihood - should this be the pseudo-likelihood instead?
    double normlihood = normLogLik(c, Ginv, littleG, curr_u, pi);

    // evaluate NB liklihood
    double c_logli = nbLogLik(mu, pc, y) - normlihood;
    double d_logli = nbLogLik(mu, pd, y) - normlihood;

    // iteratively choose either c or d and discard anything to the
    // left or right, respectively.
    while((upper - lower) >= tol){
        if(c_logli > d_logli){
            upper = pd;
            pd = pc;
            d_logli = c_logli;
            pc = lower + r*(upper - lower);
            c_logli =  nbLogLik(mu, pc, y) - normlihood;
        } else{
            lower = pc;
            pc = pd;
            c_logli = d_logli;
            pd = upper - r*(upper - lower);
            d_logli = nbLogLik(mu, pd, y) - normlihood;
        }
    }

    // return the mid-point between [c, d]
    double new_disp = (pc + pd)/2.0;
    return new_disp;
}


double phiMME(const arma::vec& y, const arma::vec& curr_sigma){
    // use the pseudovariance and method of moments
    double ps_bar = arma::var(y);
    double y_bar = arma::mean(y);
    double sigma_sum = arma::sum(y);
    double denom = 1e-2;
    double disp_mme = 1e-2;

    denom = ps_bar - y_bar - sigma_sum;
    // check for near zero denom that could blow
    // up the estimate
    if(denom < 1e-2){
        disp_mme = 1e-2;
    } else{
        disp_mme = 1/denom;
    }

    return disp_mme;
}



double nbLogLik(const arma::vec& mu, double phi, const arma::vec& y){
    // Negative binomial log-likelihood in the mean/size parameterisation, where
    // phi is the size (shape) parameter r:
    //
    //   log f(y; mu, r) = lgamma(y + r) - lgamma(r) - lgamma(y + 1)
    //                     + r * log(r / (r + mu)) + y * log(mu / (r + mu))
    //
    arma::vec denom = phi + mu;

    arma::vec logli_indiv = arma::lgamma(y + phi)
        - std::lgamma(phi)
        - arma::lgamma(y + 1.0)
        + (phi * (std::log(phi) - arma::log(denom)))
        + (y % (arma::log(mu) - arma::log(denom)));

    return arma::sum(logli_indiv);
}


double normLogLik(const int& c, const arma::mat& Ginv, const arma::mat& G,
                  const arma::vec& curr_u, double pi){
    double cdouble = (double)c;
    double detG = arma::det(G);
    double logdet = std::log(detG);

    arma::vec normlog_indiv = ((cdouble/2.0) * std::log(2*pi)) - (0.5 * logdet) - (0.5 * (curr_u.t() * Ginv * curr_u));
    double normlihood = arma::sum(normlog_indiv);

    return normlihood;
}

arma::vec heSolveREML(const arma::mat& P, const arma::vec& wdiag,
                      const Rcpp::List& dV, const arma::vec& ystar_c,
                      const bool& nnls, const bool& w_known, const int& iters){
    // Haseman-Elston estimation of the variance components from the moment
    //
    //   E[P y* y*' P] = P W P + sum_j sigma_j P dV_j P
    //
    // regressed on the vectorised lower triangles.
    //
    // How the P W P term is handled depends on whether W is known. When the
    // overdispersion is a variance component, W is the Poisson part 1/mu, fixed
    // by the current fitted means, and its coefficient is exactly 1 - so it is
    // subtracted from the response. Carrying it as a free column instead both
    // wastes information and competes with the identity basis that carries
    // sigma_0: on simulated data with a realistic spread of fitted means the two
    // correlate at 0.95, and freeing the coefficient takes the condition number
    // of the design from 5.4 to 140, which cost more than half the fits their
    // convergence.
    //
    // Under the golden section dispersion, W = (1/phi + 1/mu) carries a phi
    // estimated on a different objective. It is not known, and offsetting it
    // pushes that estimate's error straight into the components - 3 of 15 fits
    // converged and the variance component halved. There the coefficient is
    // estimated freely so the misfit is absorbed rather than propagated.
    //
    // This works off a list of partial derivatives rather than off Z and a
    // kinship matrix, so one implementation covers any set of components -
    // including the overdispersion, whose partial derivative is the identity -
    // and the same code serves fitPLGlmm, fitGeneticPLGlmm and
    // fitGeneticNullGlmm. It replaces seven near-identical estimators and four
    // vectorisers that differed only in how they assembled the same design.
    //
    // Under ML, where P is Vstar^-1 rather than the REML projection, the moment
    // carries the same mean-structure bias as any uncentred Haseman-Elston
    // estimator.
    const int ctot = dV.size();
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(P));

    arma::vec Py = P * ystar_c;
    arma::mat Ycov = Py * Py.t();

    // W is diagonal by construction, so keep the first product O(n^2)
    arma::mat wbasis = (P.each_row() % wdiag.t()) * P.t();

    const int ncol = w_known ? ctot : ctot + 1;
    const int off = w_known ? 0 : 1;

    arma::vec Ybig = Ycov(lower_indices);
    if(w_known){
        Ybig -= wbasis(lower_indices);
    }

    arma::mat design(Ybig.n_elem, ncol);
    if(!w_known){
        design.col(0) = wbasis(lower_indices);
    }
    for(int j = 0; j < ctot; j++){
        const arma::mat& dVj = dV[j];
        arma::mat basis = (P * dVj) * P.t();
        design.col(j + off) = basis(lower_indices);
    }

    arma::vec est(ncol, arma::fill::zeros);
    if(nnls){
        est = nnlsSolve(design, Ybig, est, iters);
    } else{
        est = arma::solve(design, Ybig, arma::solve_opts::fast);
    }

    return est.tail(ctot);
}
