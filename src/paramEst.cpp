#include "paramEst.h"
#include "computeMatrices.h"
#include "utils.h"
#include<RcppArmadillo.h>
#include<RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// using namespace Rcpp;

// All functions used in parameter estimation

arma::vec sigmaScoreREML_arma (Rcpp::List pvstar_i, const arma::vec& ystar, const arma::mat& P){
    // Armadillo implementation
    const int& c = pvstar_i.size();
    arma::vec reml_score(c);
    // arma::sp_mat _P(P);

    for(int i=0; i < c; i++){
        const arma::mat& P_pvi = pvstar_i(i); // this is P * partial derivative

        double lhs = -0.5 * arma::trace(P_pvi);
        arma::mat mid1(1, 1);
        mid1 = arma::trans(ystar) * P_pvi * P * ystar;
        double rhs = 0.5 * mid1[0, 0];

        reml_score[i] = lhs + rhs;
    }

    return reml_score;
}


arma::mat sigmaInfoREML_arma (const Rcpp::List& pvstari, const arma::mat& P){
    // REML Fisher/expected information matrix
    const int& c = pvstari.size();
    arma::mat sinfo(c, c);

    // this is a symmetric matrix so only need to fill the upper or
    // lower triangle to make it O(n^2/2) rather than O(n^2)
    for(int i=0; i < c; i++){
        const arma::mat& _ipP = pvstari[i]; // this is P * \d Var/ \dsigma

        for(int j=i; j < c; j++){
            const arma::mat& P_jp = pvstari[j]; // this is P * \d Var/ \dsigma
            arma::mat a_ij(_ipP * P_jp); // this is the biggest bottleneck - it takes >2s!
            double _artr = arma::trace(a_ij);

            sinfo(i, j) = 0.5 * _artr;
            if(i != j){
                sinfo(j, i) = 0.5 * _artr;
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


arma::vec fisherScore (arma::mat hess, arma::vec score_vec, arma::vec theta_hat){
    // sequentially update the parameter using the Newton-Raphson algorithm
    // theta ~= theta_hat + hess^-1 * score
    // this needs to be in a direction of descent towards a minimum
    int m = theta_hat.size();
    arma::vec theta(m);
    arma::mat hessinv(hess.n_rows, hess.n_cols);

    // will need a check here for singular hessians...
    try{
        double _rcond = arma::rcond(hess);
        bool is_singular;
        is_singular = _rcond < 1e-9;

        // check for singular condition
        if(is_singular){
            Rcpp::stop("Variance Component Hessian is computationally singular");
        }

        hessinv = arma::inv(hess); // always use pinv? solve() and inv() are most sensitive than R versions
        theta = theta_hat + (hessinv * score_vec);
        return theta;
    } catch(std::exception &ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }
}


arma::mat coeffMatrix(const arma::mat& X, const arma::mat& Winv, const arma::mat& Z, const arma::mat& Ginv){
    // compute the components of the coefficient matrix for the MMEs
    //
    int c = Z.n_cols;
    int m = X.n_cols;

    arma::mat ul(m, m);
    arma::mat ur(m, c);
    arma::mat ll(c, m);
    arma::mat lr(c, c);

    arma::mat lhs_top(m, m+c);
    arma::mat lhs_bot(c, m+c);
    arma::mat lhs(m+c, m+c);

    ul = X.t() * Winv * X;
    ur = X.t() * Winv * Z;
    ll = Z.t() * Winv * X;
    lr = (Z.t() * Winv * Z) + Ginv;

    lhs_top = arma::join_rows(ul, ur); // join_rows matches the rows i.e. glue columns together
    lhs_bot = arma::join_rows(ll, lr);

    lhs = arma::join_cols(lhs_top, lhs_bot); // join_cols matches the cols, i.e. glue rows together
    return lhs;
}


arma::vec solveEquations (const int& c, const int& m, const arma::mat& Winv, const arma::mat& Zt, const arma::mat& Xt,
                          arma::mat coeffmat, arma::vec beta, arma::vec u, const arma::vec& ystar){
    // solve the mixed model equations
    arma::vec rhs_beta(m);
    arma::vec rhs_u(c);
    arma::mat rhs(m+c, 1);

    arma::vec theta_up(m+c);

    rhs_beta.col(0) = Xt * Winv * ystar;
    rhs_u.col(0) = Zt * Winv * ystar;

    rhs = arma::join_cols(rhs_beta, rhs_u);

    // need a check for singular hessian here
    try{
        double _rcond = arma::rcond(coeffmat);
        bool is_singular;
        is_singular = _rcond < 1e-9;

        // check for singular condition
        if(is_singular){
            // this happens when G^-1 contains NaN values <- how does this happen
            // and how do we prevent it?
            Rcpp::Rcout << _rcond << std::endl;
            // Rcpp::stop("Coefficients Hessian is computationally singular");
            throw std::range_error("Coefficients Hessian is computationally singular");
        }

        // can we just use solve here instead?
        // if the coefficient matrix is singular then do we resort to pinv?
        theta_up = arma::solve(coeffmat, rhs);

        return theta_up;
    } catch(std::exception &ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }

}


arma::mat computeZstar(const arma::mat& Z, arma::vec curr_sigma, Rcpp::List u_indices){
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
        curr_sigma.brief_print("Sigma\n");
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


arma::vec solveEquationsPCG (const int& c, const int& m, const arma::mat& Winv, const arma::mat& Zt, const arma::mat& Xt,
                          arma::mat coeffmat, arma::vec curr_theta, const arma::vec& ystar, const double& conv_tol){
    // solve the mixed model equations with a preconditioned conjugate gradient
    // A = coeffmat
    // x = theta
    // b = rhs

    arma::vec rhs_beta(m);
    arma::vec rhs_u(c);
    arma::mat rhs(m+c, 1);

    arma::vec theta_up(m+c);

    rhs_beta.col(0) = Xt * Winv * ystar;
    rhs_u.col(0) = Zt * Winv * ystar;

    rhs = arma::join_cols(rhs_beta, rhs_u);

    // I'll assume any preconditioning has already been applied
    // need a check for singular hessian here
    try{
        double _rcond = arma::rcond(coeffmat);
        bool is_singular;
        is_singular = _rcond < 1e-9;

        // check for singular condition
        if(is_singular){
            Rcpp::stop("Coefficients Hessian is computationally singular");
        }

        // can we just use solve here instead?
        // if the coefficient matrix is singular then do we resort to pinv?
        // is it worth doing a quick analysis of the eigen values of coeff?
        // the _rcond might be sufficient to tell us if the matrix is ill-conditioned
        // do we need to know a priori if we have a few large eigenvalues??
        // maybe this could be tweaked by setting the convergence criterai to > 0?
        theta_up = conjugateGradient(coeffmat, curr_theta, rhs, conv_tol);
        // theta_up = arma::solve(coeffmat, rhs);

        return theta_up;
    } catch(std::exception &ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }

}


arma::vec conjugateGradient(arma::mat A, arma::vec x, arma::vec b, double conv_tol){
    // use conjugate gradients to solve the system of linear equations
    // Ax = b
    // Algorithm:
    // r_0 <- Ax_0 - b, p_0 <- -r_0, k <- 0
    // while r_k != 0:
    // alpha_k <- (rk^T * rk)/(pk^T * A * pk)
    // x_k+1 <- xk + alpha_k * pK
    // r_k+1 <- rk + alpha_k * A * pK
    // beta_k+1 <- r rk+1 + beta_k+1 * pk
    // pk+1 <- -r_k+1 + beta_k+1 * pk
    // k++

    // need to generate x_0 from the current estimates: [beta u]
    const unsigned int m = A.n_cols;
    const unsigned int n = b.size();

    arma::dcolvec xk(m);
    xk = x; // use current estimates as x0
    arma::vec xk_update(m);
    xk_update = arma::dcolvec(m);
    // x0.randu(); // initial x values

    double alpha_k = 0.0;
    double beta_k = 0.0;

    arma::dcolvec rk(n);
    arma::dcolvec rk_update(n);
    arma::dcolvec pk(m);
    arma::dcolvec pk_update(m);

    rk = (A * xk) - b;
    pk = -rk;
    unsigned int k = 0;

    Rcpp::LogicalVector _check_rzero = check_tol_arma_numeric(rk, conv_tol);
    bool _all_rk_zero = Rcpp::all(_check_rzero).is_false(); // .is_false() required for proper type casting to bool

    while(_all_rk_zero){ // evaluates true until all rk are zero
        alpha_k = (rk.t() * rk).eval()(0,0)/(pk.t() * A * pk).eval()(0, 0); // needed to convert vector inner product to scalar
        xk_update = xk + alpha_k * pk;
        rk_update = rk + alpha_k * A * pk;
        beta_k = (rk_update.t() * rk_update).eval()(0, 0)/(rk.t() * rk).eval()(0, 0); // needed to convert vector inner product to scalar
        pk_update = -rk_update + beta_k * pk;

        rk = rk_update;
        pk = pk_update;
        xk = xk_update;

        _check_rzero = check_tol_arma_numeric(rk, conv_tol);
        _all_rk_zero = Rcpp::all(_check_rzero).is_false(); // .is_false() required for proper type casting to bool
        k++;
    }

    Rprintf("CG completed in %u iterations", k);
    return xk_update;
}
