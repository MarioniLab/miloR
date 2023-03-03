#include "paramEst.h"
#include "computeMatrices.h"
#include "utils.h"
#include<RcppArmadillo.h>
#include<RcppEigen.h>
#include<cmath>
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
    arma::vec theta(m, arma::fill::zeros);
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
    } catch(std::exception const& ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }

    return theta;
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

    arma::vec theta_up(m+c, arma::fill::zeros);

    rhs_beta.col(0) = Xt * Winv * ystar;
    rhs_u.col(0) = Zt * Winv * ystar;

    rhs = arma::join_cols(rhs_beta, rhs_u);

    // need a check for singular hessian here
    double _rcond = arma::rcond(coeffmat);
    try{
        // double _rcond = arma::rcond(coeffmat);
        bool is_singular;
        is_singular = _rcond < 1e-9;

        // check for singular condition
        if(is_singular){
            // this happens when G^-1 contains NaN values <- how does this happen
            // and how do we prevent it?
            Rcpp::Rcout << "Condition number: " << _rcond << std::endl;
            throw std::runtime_error("Coefficients Hessian is computationally singular");
        }

        // can we just use solve here instead?
        // if the coefficient matrix is singular then do we resort to pinv?
        theta_up = arma::solve(coeffmat, rhs, arma::solve_opts::no_approx);
    } catch(std::exception const& ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }
    return theta_up;
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

    arma::vec theta_up(m+c, arma::fill::zeros);

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

    } catch(std::exception &ex){
        forward_exception_to_r(ex);
    } catch(...){
        Rf_error("c++ exception (unknown reason)");
    }

    return theta_up;
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

    Rprintf("CG completed in %u iterations\n", k);
    return xk_update;
}


arma::vec estHasemanElstonGenetic(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices, arma::vec ystar, arma::mat Kin){
    // use HasemanElston regression to estimate variance components
    // vectorize everything
    // we will also estimate a "residual" variance parameter
    unsigned int n = ystar.size();
    unsigned int c = u_indices.size(); // number of variance components
    unsigned long nsq = n * (n + 1)/2; //size of vectorised components using just upper or lower triangle of covariance matrix

    arma::mat Ycovar(n, n);
    Ycovar = PREML * (ystar * ystar.t()) * PREML; // project onto REML P matrix

    // select the upper triangular elements, including the diagonal
    arma::uvec upper_indices = trimatu_ind(arma::size(Ycovar));
    arma::vec Ybig = Ycovar(upper_indices);

    // sequentially vectorise ZZ^T - this automatically adds a vectorised identity matrix
    // for the "residual" variance
    arma::mat vecZ(nsq, c+1);
    vecZ = vectoriseZGenetic(Z, u_indices, PREML, Kin); // projection already applied

    // solve by linear least squares
    arma::vec _he_update(c+1);
    arma::vec he_update(c);

    // use OSL here for starters
    _he_update = arma::solve(vecZ, Ybig);
    he_update = _he_update.tail(c);

    return he_update;
}


arma::vec estHasemanElston(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices, arma::vec ystar){
    // use HasemanElston regression to estimate variance components
    // vectorize everything
    // we will also estimate a "residual" variance parameter
    unsigned int n = ystar.size();
    unsigned int c = u_indices.size(); // number of variance components
    unsigned long nsq = n * (n + 1)/2; //size of vectorised components using just upper or lower triangle of covariance matrix

    arma::mat Ycovar(n, n);
    Ycovar = PREML * (ystar * ystar.t()) * PREML; // project onto REML P matrix

    // select the upper triangular elements, including the diagonal
    arma::uvec upper_indices = trimatu_ind(arma::size(Ycovar));
    arma::vec Ybig = Ycovar(upper_indices);

    // sequentially vectorise ZZ^T - this automatically adds a vectorised identity matrix
    // for the "residual" variance
    arma::mat vecZ(nsq, c+1);
    vecZ = vectoriseZ(Z, u_indices, PREML); // projection already applied

    // solve by linear least squares
    arma::vec _he_update(c+1);
    arma::vec he_update(c);

    // use OSL here for starters
    _he_update = arma::solve(vecZ, Ybig);
    he_update = _he_update.tail(c);

    return he_update;
}


arma::vec estHasemanElstonConstrained(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices,
                                      arma::vec ystar, arma::vec he_update, const int& Iters){
    // use constrained HasemanElston regression to estimate variance components - using a NNLS estimator
    // vectorize everything
    // we will also estimate a "residual" variance parameter
    // however, there is no reason this "residual" paramer has to be constrained...
    unsigned int n = ystar.size();
    unsigned int c = u_indices.size(); // number of variance components
    unsigned long nsq = n * (n + 1)/2; //size of vectorised components using just upper or lower triangle of covariance matrix

    arma::mat Ycovar(n, n);
    Ycovar = PREML * (ystar * ystar.t()) * PREML; // project onto REML P matrix

    // select the upper triangular elements, including the diagonal
    arma::uvec upper_indices = trimatu_ind(arma::size(Ycovar));
    arma::vec Ybig = Ycovar(upper_indices);

    // sequentially vectorise ZZ^T - this automatically adds a vectorised identity matrix
    // for the "residual" variance
    arma::mat vecZ(nsq, c+1);
    vecZ = vectoriseZ(Z, u_indices, PREML); // projection already applied

    // solve by linear least squares
    arma::vec _he_update(c+1);

    // first check if we can get a non-negative OLS estimate
    arma::dvec _ols = arma::solve(vecZ, Ybig, arma::solve_opts::fast);
    if(any(_ols < 1e-8)){
        // use NNSL here - Lawson and Hanson algorithm or FAST-NNLS
        // the latter is only applicable when vecZ is PD - need to check this with all positive eigenvalues
        bool _ispd;
        _ispd = check_pd_matrix(vecZ);

        if(_ispd){
            // use the FAST NNLS solver
            _he_update = fastNnlsSolve(vecZ, Ybig);
        } else{
            // have to use slower implementation from Lawson and Hanson
            _he_update = nnlsSolve(vecZ, Ybig, he_update, Iters);
        }
        he_update = _he_update.tail(c);
    } else{
        he_update = _ols.tail(c);
    }

    return he_update;
}


arma::vec estHasemanElstonConstrainedGenetic(const arma::mat& Z, const arma::mat& PREML, Rcpp::List u_indices,
                                             arma::vec ystar, arma::mat Kin, arma::vec he_update, const int& Iters){
    // use constrained HasemanElston regression to estimate variance components - using a NNLS estimator
    // vectorize everything
    // we will also estimate a "residual" variance parameter
    unsigned int n = ystar.size();
    unsigned int c = u_indices.size(); // number of variance components
    unsigned long nsq = n * (n + 1)/2; //size of vectorised components using just upper or lower triangle of covariance matrix

    arma::mat Ycovar(n, n);
    Ycovar = PREML * (ystar * ystar.t()) * PREML; // project onto REML P matrix

    // select the upper triangular elements, including the diagonal
    arma::uvec upper_indices = trimatu_ind(arma::size(Ycovar));
    arma::vec Ybig = Ycovar(upper_indices);

    // sequentially vectorise ZZ^T - this automatically adds a vectorised identity matrix
    // for the "residual" variance
    arma::mat vecZ(nsq, c+1);
    vecZ = vectoriseZGenetic(Z, u_indices, PREML, Kin); // projection already applied

    arma::vec _he_update(c+1);

    // first check if we can get a non-negative OLS estimate
    arma::dvec _ols = arma::solve(vecZ, Ybig, arma::solve_opts::fast);
    if(any(_ols < 1e-8)){
        // use NNSL here - Lawson and Hanson algorithm or FAST-NNLS
        // the latter is only applicable when vecZ is PD - need to check this with all positive eigenvalues
        bool _ispd;
        _ispd = check_pd_matrix(vecZ);

        if(_ispd){
            // use the FAST NNLS solver
            _he_update = fastNnlsSolve(vecZ, Ybig);
        } else{
            // have to use slower implementation from Lawson and Hanson
            _he_update = nnlsSolve(vecZ, Ybig, he_update, Iters);
        }
        he_update = _he_update.tail(c);
    } else{
        he_update = _ols.tail(c);
    }

    return he_update;
}


arma::vec nnlsSolve(const arma::mat& vecZ, arma::vec Y, arma::vec nnls_update, const int& Iters){
    // Lawson and Hanson algorithm for constrained NNLS
    // initial params
    // P is the passive set - i.e. the indices of the non-negative estimates
    // R is the active set - the indices of the constrained estimates, i.e. those held at zero
    // d = 0 is a solution vector
    // w is a vector of the Langrangian multipliers

    // need to implement a check for infinite loop - have a max iters argument?
    double constval = 0.0; // value at which to constrain values
    unsigned int inner_count = 0;
    unsigned int outer_count = 0;
    double EPS = 1e-6;
    unsigned int m = vecZ.n_cols;
    arma::ivec P(m, arma::fill::ones); // the indices have to be set to negative values to be empty

    for(int i=0; i < m; i++){
        if(nnls_update[i] <= constval){
            P[i] = -1;
        }
    }

    // all indices need to be active if all are zero - i.e. the first iteration
    // what happens if all estimates get regularized to exactly zero?
    arma::ivec R(m, arma::fill::ones); // these are the indices of the active set

    for(int i=0; i < m; i++){
        if(nnls_update[i] > constval){
            R[i] = -1;
        }
    }

    arma::dvec w(m);
    w = vecZ.t() * (Y - vecZ*nnls_update); //Lagrangian multipliers

    // arma::mat _hessian = arma::inv(vecZ.t() * vecZ);
    // _hessian.print("(X^TX)^-1");
    bool check_R;
    bool check_wR;
    bool check_wP;
    bool check_conditions;
    check_R = all(R < 0);
    check_wR = all(w.elem(find(R > 0)) <= constval); // these should be 0 or negative
    check_wP = all(abs(w.elem(find(P > 0)) - constval) < EPS); // these should be ~0
    check_conditions = check_R || (check_wR && check_wP);

    double max_w = 0.0;

    while(!check_conditions){ // set a tolerance here to check for positive Langrangian
        max_w = max(w.elem(find(R > 0))); // what if several are identical?
        unsigned int max_j = 0;

        // get the index of the maximal Langrangian multiplier
        // need to find the maximum value in w^R, but the index needs to be from w
        // turn this into a function?
        for(int i=0; i < m; i++){
            if(abs(w[i] - max_w) <= EPS){
                max_j = i;
            }
        }

        P[max_j] = -1 * P[max_j]; // this should reverse the sign
        R[max_j] = -1 * R[max_j];

        // find the elements >= 0
        arma::uvec select_P = find(P > 0);
        arma::vec s_all(m);
        s_all.fill(constval);
        s_all.elem(select_P) = arma::solve(vecZ.cols(select_P), Y); // is this faster?

        double min_sp = s_all.elem(select_P).min();
        double alpha; // the step size
        // outer_count++;

        while(min_sp <= 0){
            inner_count = 0;
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
            s_all.elem(select_P) = arma::solve(vecZ.cols(select_P), Y);
            min_sp = s_all.elem(select_P).min();
            // inner_count++;
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


arma::vec fastNnlsSolve(const arma::mat& vecZ, arma::vec Y){
    // This uses the Fast Approximate Solution Trajectory NNLS
    // from https://stackoverflow.com/questions/58006606/rcpp-implementation-of-fast-nonnegative-least-squares

    // assumes the problem is not ill-conditioned - how realistic is this?
    // vecZ is almost _never_ square
    int m = vecZ.n_rows;
    arma::mat A(m, m);
    A = vecZ * vecZ.t(); // is this right? I don't think that it is!
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


arma::mat vectoriseZ(arma::mat Z, Rcpp::List u_indices, arma::mat P){
    // sequentially vectorise the columns ZZ^T that map to each random effect
    // pre- and post- multiply by the REML projection matrix
    int c = u_indices.size();
    int n = Z.n_rows;
    unsigned long nsq = n * (n + 1)/2;

    Rcpp::List _Zelements(1);
    arma::mat bigI(n, n, arma::fill::eye); // this should be the identity which we vectorise

    // select the upper triangular elements, including the diagonal
    arma::uvec upper_indices = arma::trimatu_ind(arma::size(bigI));
    arma::mat _vI = bigI(upper_indices);
    _Zelements(0) = _vI;

    for(int i=0; i < c; i++){
        // extract the elements of u_indices
        arma::uvec u_idx = u_indices(i);
        unsigned int q = u_idx.n_rows;
        arma::mat _subZ(n, q);
        unsigned int qmin = arma::min(u_idx);
        unsigned int qmax = arma::max(u_idx);
        _subZ = Z.cols(qmin-1, qmax-1);

        arma::mat _ZZT(n, n);
        _ZZT = P * (_subZ * _subZ.t()) * P; // REML projection

        // vectorise
        arma::vec _vecZ = _ZZT(upper_indices);
        arma::mat _vecZZT = _Zelements(0);
        unsigned long _zc = _vecZZT.n_cols;

        arma::mat _vZ(nsq, _zc+1);
        _vZ = arma::join_rows(_vecZZT, _vecZ);
        _Zelements(0) = _vZ;
    }

    arma::mat vecMat = _Zelements(0);
    return vecMat;
}


arma::mat vectoriseZGenetic(arma::mat Z, Rcpp::List u_indices, arma::mat P, arma::mat Kin){
    // sequentially vectorise the columns ZZ^T that map to each random effect
    // pre- and post- multiply by the REML projection matrix
    int c = u_indices.size();
    int n = Z.n_rows;
    // unsigned long nsq = pow(static_cast<float>(n), 2);
    unsigned long nsq = n * (n + 1)/2;

    Rcpp::List _Zelements(1);
    arma::mat bigI(n, n, arma::fill::ones); // should this be column of 1s

    // select the upper triangular elements, including the diagonal
    arma::uvec upper_indices = arma::trimatu_ind(arma::size(bigI));
    arma::mat _vI = bigI(upper_indices);
    _Zelements(0) = _vI;

    for(int i=0; i < c; i++){
        // extract the elements of u_indices
        arma::uvec u_idx = u_indices(i);
        unsigned int q = u_idx.n_rows;
        arma::mat _subZ(n, q);
        unsigned int qmin = arma::min(u_idx);
        unsigned int qmax = arma::max(u_idx);
        _subZ = Z.cols(qmin-1, qmax-1);

        // always set the last component to the genetic variance if there is a kinship matrix
        if(i == c-1){
            arma::vec _vecZ = Kin(upper_indices);
            arma::mat _vecZZT = _Zelements(0);
            unsigned long _zc = _vecZZT.n_cols;

            arma::mat _vZ(nsq, _zc+1);
            _vZ = arma::join_rows(_vecZZT, _vecZ);
            _Zelements(0) = _vZ;
        } else{
            // compute Z_i Z_i^T
            arma::mat _ZZT(n, n);
            _ZZT = P * (_subZ * _subZ.t()) * P; // REML projection

            // vectorise
            arma::vec _vecZ = _ZZT(upper_indices);
            arma::mat _vecZZT = _Zelements(0);
            unsigned long _zc = _vecZZT.n_cols;

            arma::mat _vZ(nsq, _zc+1);
            _vZ = arma::join_rows(_vecZZT, _vecZ);
            _Zelements(0) = _vZ;
        }
    }

    arma::mat vecMat = _Zelements(0);
    return vecMat;
}


double phiLineSearch(double disp, double lower, double upper, const int& c,
                     arma::vec mu, arma::mat Ginv, double pi,
                     arma::vec curr_u, arma::vec sigma, arma::vec y){
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
                       arma::vec mu, arma::mat Ginv, double pi,
                       arma::vec curr_u, arma::vec sigma, arma::vec y){

    // perform a golden section search
    // function inspired by https://drlvk.github.io/nm/section-golden-section.html
    // this only requires a single call to this function

    double tol = 1e-2; // tolerance for phi
    double r = (3 - std::sqrtf(5))/2.0;
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



double nbLogLik(arma::vec mu, double phi, arma::vec y){
    double logli = 0.0;
    arma::vec logli_indiv(y.n_rows);
    arma::vec muphi(y.n_rows);
    muphi = mu/(mu + phi);

    // element wise multiplication of y and other equation elements
    logli_indiv = y % arma::log(mu/(mu + phi)) + (phi * (1 - (mu/(mu+phi)))) + (arma::lgamma(y+1) - std::lgamma(phi));

    logli = arma::sum(logli_indiv);
    return logli;
}


double normLogLik(const int& c, arma::mat Ginv, arma::mat G, arma::vec curr_u, double pi){
    double cdouble = (double)c;
    double detG = arma::det(G);
    double logdet = std::log(detG);

    arma::vec normlog_indiv = ((cdouble/2.0) * std::log(2*pi)) - (0.5 * logdet) - (0.5 * (curr_u.t() * Ginv * curr_u));
    double normlihood = arma::sum(normlog_indiv);

    return normlihood;
}



