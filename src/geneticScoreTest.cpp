#include<RcppArmadillo.h>
#include<string>
#include<cmath>
#include<limits>
// [[Rcpp::depends(RcppArmadillo)]]
#include "paramEst.h"
#include "computeMatrices.h"
#include "utils.h"
#include "geneticScoreTest.h"
using namespace Rcpp;


// Build the pseudo-variance V* = diag(w) + sum_j sigma_j Z_j Z_j' + sigma_c K
//
// This never materialises a dense n x n block inside Z for the genetic random
// effect. The genetic term enters directly as sigma_c * K, which is the
// Z_g = I, G_g = sigma_c * K parameterisation. That is algebraically identical
// to Z_g = chol(K), G_g = sigma_c * I but avoids an n x n dense Z block, and it
// applies K exactly once.
arma::mat buildVstar(const arma::vec& wdiag, const arma::mat& Z, const arma::mat& K,
                     const arma::vec& sigmas, const Rcpp::List& u_indices){
    const unsigned int n = wdiag.n_elem;
    const int c = sigmas.n_elem;

    arma::mat Vstar(n, n, arma::fill::zeros);
    Vstar.diag() = wdiag;

    // non-genetic random effects: sigma_j * Z_j Z_j'
    for(int j = 0; j < c - 1; j++){
        arma::uvec jdx = u_indices[j];
        arma::mat Zj = Z.cols(jdx - 1); // R is 1-based
        Vstar += sigmas(j) * (Zj * Zj.t());
    }

    // genetic random effect
    Vstar += sigmas(c - 1) * K;

    // enforce exact symmetry - guards inv_sympd against accumulated round-off
    Vstar = 0.5 * (Vstar + Vstar.t());
    return Vstar;
}


// Partial derivatives of V* wrt each variance component. For the non-genetic
// effects this is Z_j Z_j'; for the genetic effect it is K. These do not depend
// on the variance components, so they are built once per fit rather than per
// iteration.
Rcpp::List buildVpartial(const arma::mat& Z, const arma::mat& K, const Rcpp::List& u_indices,
                         const int& c){
    Rcpp::List dV(c);
    for(int j = 0; j < c - 1; j++){
        arma::uvec jdx = u_indices[j];
        arma::mat Zj = Z.cols(jdx - 1);
        dV[j] = arma::mat(Zj * Zj.t());
    }
    dV[c - 1] = K;
    return dV;
}


//' Fit a NB-GLMM with a genetic random effect, optionally warm-started
//'
//' Fits the negative binomial generalised linear mixed model used by Milo for
//' genetic analyses, using penalised quasi-likelihood (PQL) with REML or ML
//' estimation of the variance components. It follows the same iterative logic as
//' \code{fitGeneticPLGlmm} - alternating between the fixed/random effect
//' solutions and the variance component update, with an off-line golden-section
//' search for the dispersion - but differs in three respects that matter for
//' genome-wide use:
//'
//' \enumerate{
//' \item The genetic covariance enters as \eqn{\sigma_g K} directly rather than
//' through a dense \eqn{n \times n} Cholesky block appended to \emph{Z}. This
//' applies \emph{K} exactly once and removes the dominant \eqn{O(n^3)} term.
//' \item The working weight matrices \emph{D} and \emph{W} are diagonal by
//' construction and are carried as vectors, so they are never inverted by a
//' dense LU factorisation.
//' \item Null model parameter estimates may be supplied, in which case they are
//' used as starting values and - if \code{fix_variance} is \code{TRUE} - the
//' variance components and dispersion are held fixed while only the fixed
//' effects are updated.
//' }
//'
//' When \code{return_projection} is \code{TRUE} the REML projection matrix
//' \emph{P} and the vector \eqn{P y^*} are returned. These are the sufficient
//' quantities for the score test in \code{\link{scoreTestGeneticSNPs}}, and are
//' invariant across the SNPs tested within a neighbourhood.
//'
//' @param Z mat - design matrix mapping the levels of the non-genetic random
//' effects to observations. May have zero columns if the genetic random effect
//' is the only one.
//' @param X mat - design matrix of fixed effects. For a null model this excludes
//' the genetic variant of interest.
//' @param K mat - the n x n genetic relationship matrix.
//' @param muvec vec - vector of starting values for the phenotype means.
//' @param offsets vec - vector of model offsets, e.g. log library sizes.
//' @param curr_beta vec - starting values for the fixed effect parameters.
//' @param curr_u vec - starting values for the random effect BLUPs, ordered as
//' the non-genetic effect levels followed by the n genetic effect levels.
//' @param curr_sigma vec - starting values for the variance components. The
//' final element is always the genetic variance component attached to \emph{K}.
//' @param y vec - vector of observed counts.
//' @param u_indices List - each element holds the (1-based) column indices of
//' \emph{Z} belonging to one non-genetic random effect.
//' @param theta_conv double - convergence tolerance on the parameter estimates.
//' @param curr_disp double - starting value for the dispersion.
//' @param REML bool - use REML rather than ML for the variance components.
//' @param maxit int - maximum number of PQL iterations.
//' @param Kinv_ Nullable NumericMatrix - optional precomputed inverse of
//' \emph{K}. \emph{K} is invariant across neighbourhoods and SNPs, so supplying
//' the inverse once avoids repeating an O(n^3) factorisation per model fit.
//' @param null_sigma_ Nullable NumericVector - optional variance components from
//' a previously fitted null model, used as starting values.
//' @param null_beta_ Nullable NumericVector - optional fixed effect estimates
//' from a previously fitted null model, used as starting values. Only the
//' leading elements shared with the current design are used.
//' @param null_disp double - optional dispersion estimate from a previously
//' fitted null model. A negative value means no estimate is supplied.
//' @param fix_variance bool - hold the variance components and dispersion fixed
//' at their supplied values and update only the fixed effects.
//' @param return_projection bool - return the REML projection matrix \emph{P}
//' and \eqn{P y^*}.
//' @param fix_dispersion bool - hold the dispersion fixed at \code{curr_disp}
//' rather than re-estimating it inside the PQL loop. Recommended when a
//' dispersion estimated across neighbourhoods is available, e.g. from
//' \code{edgeR::estimateDisp}.
//' @param max_disp double - upper bound on the size parameter when it is
//' estimated inside the loop. Prevents the search running away to the Poisson
//' limit.
//' @param disp_as_vc bool - estimate the negative binomial overdispersion as an
//' additional variance component on the REML objective, rather than by a
//' golden-section search on the conditional negative binomial likelihood. This
//' puts the overdispersion and the random effect variances on a common
//' objective so that they compete properly. Defaults to \code{TRUE}; set to
//' \code{FALSE} to use the golden-section search on the conditional negative
//' binomial likelihood instead.
//'
//' @details The model fitted is the same pseudo-likelihood approximation used
//' throughout Milo. At convergence the working response is
//' \eqn{y^* = \eta + D^{-1}(y - \mu)} and the pseudo-variance is
//' \eqn{V^* = W + \sum_j \sigma_j Z_j Z_j' + \sigma_g K}, with
//' \eqn{W = \mathrm{diag}(\phi^{-1} + \mu_i^{-1})}. Variance components are
//' updated by Fisher scoring on the REML log-likelihood and constrained to be
//' non-negative.
//'
//' @return A \code{list} containing the fitted model. Elements are:
//' \describe{
//' \item{\code{FE}:}{\code{numeric} vector of fixed effect estimates.}
//' \item{\code{RE}:}{\code{numeric} vector of random effect BLUPs.}
//' \item{\code{Sigma}:}{\code{numeric} vector of variance component estimates, genetic component last.}
//' \item{\code{Dispersion}:}{\code{numeric} scalar dispersion estimate.}
//' \item{\code{converged}:}{\code{logical} whether the convergence tolerance was met.}
//' \item{\code{Iters}:}{\code{numeric} number of PQL iterations run.}
//' \item{\code{SE}:}{\code{numeric} vector of fixed effect standard errors.}
//' \item{\code{t}:}{\code{numeric} vector of Wald t-statistics for the fixed effects.}
//' \item{\code{P}:}{\code{matrix} REML projection matrix, or a 1 x 1 zero matrix if not requested.}
//' \item{\code{Pystar}:}{\code{numeric} vector \eqn{P y^*}, or a length-1 zero vector if not requested.}
//' \item{\code{ystar}:}{\code{numeric} working response at convergence.}
//' \item{\code{Wdiag}:}{\code{numeric} diagonal of the working weight matrix.}
//' \item{\code{VCOV}:}{\code{matrix} variance-covariance matrix of the fixed effects.}
//' \item{\code{LOGLIHOOD}:}{\code{numeric} pseudo-log-likelihood at convergence.}
//' }
//'
//' @author Mike Morgan
//'
//' @examples
//' NULL
//'
//' @name fitGeneticNullGlmm
//'
// [[Rcpp::export]]
List fitGeneticNullGlmm(const arma::mat& Z, const arma::mat& X, const arma::mat& K,
                        arma::vec muvec, arma::vec offsets,
                        arma::vec curr_beta, arma::vec curr_u, arma::vec curr_sigma,
                        const arma::vec& y, List u_indices,
                        double theta_conv, double curr_disp,
                        const bool& REML, const int& maxit,
                        Rcpp::Nullable<Rcpp::NumericMatrix> Kinv_ = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> null_sigma_ = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> null_beta_ = R_NilValue,
                        double null_disp = -1.0,
                        const bool& fix_variance = false,
                        const bool& return_projection = true,
                        const bool& fix_dispersion = false,
                        double max_disp = 1e4,
                        const bool& disp_as_vc = true){

    constexpr double pi = 3.14159265358979323846;
    const double constval = 1e-8;

    const int c = curr_sigma.n_elem;
    const int m = X.n_cols;
    const unsigned int n = X.n_rows;
    const int qz = Z.n_cols;

    if(K.n_rows != n || K.n_cols != n){
        stop("Genetic relationship matrix dimensions do not match the number of observations");
    }
    if(c < 1){
        stop("At least one variance component is required");
    }

    // ---- optional warm start from a previously fitted null model -----------
    if(null_sigma_.isNotNull()){
        arma::vec _ns = as<arma::vec>(NumericVector(null_sigma_));
        if(_ns.n_elem == static_cast<unsigned int>(c)){
            curr_sigma = _ns;
        }
    }

    if(null_beta_.isNotNull()){
        arma::vec _nb = as<arma::vec>(NumericVector(null_beta_));
        // the alternative design shares its leading columns with the null design
        unsigned int nshare = std::min(_nb.n_elem, curr_beta.n_elem);
        for(unsigned int i = 0; i < nshare; i++){
            curr_beta(i) = _nb(i);
        }
    }

    if(null_disp > 0.0){
        curr_disp = null_disp;
    }

    // ---- K inverse: invariant across neighbourhoods and SNPs ---------------
    arma::mat Kinv(n, n);
    if(Kinv_.isNotNull()){
        Kinv = as<arma::mat>(NumericMatrix(Kinv_));
    } else {
        bool _kok = arma::inv_sympd(Kinv, K);
        if(!_kok){
            Rcpp::warning("Genetic relationship matrix is not positive definite - using pseudoinverse");
            Kinv = arma::pinv(K);
        }
    }

    // Partial derivatives of V* are constant across iterations.
    //
    // When disp_as_vc is set, the negative binomial overdispersion is treated as
    // one more variance component. W = diag(1/phi + 1/mu) splits into a constant
    // diagonal 1/phi and the Poisson part 1/mu, and the constant diagonal is
    // exactly sigma_0 * I. Estimating sigma_0 = 1/phi by the same REML Fisher
    // scoring as the other components puts every variance parameter on one
    // objective, instead of estimating phi from the conditional NB likelihood at
    // fixed mu while sigma is estimated from the marginal REML likelihood.
    const int ctot = disp_as_vc ? c + 1 : c;
    List dV = buildVpartial(Z, K, u_indices, c);
    List dVa(ctot);
    for(int j = 0; j < c; j++){
        dVa[j] = dV[j];
    }
    if(disp_as_vc){
        dVa[c] = arma::mat(n, n, arma::fill::eye);
    }

    // augmented parameter vector: [sigma_1 .. sigma_c, sigma_0]
    arma::vec sig_a(ctot);
    sig_a.head(c) = curr_sigma;
    if(disp_as_vc){
        sig_a(c) = 1.0 / std::max(curr_disp, 1e-8);
    }

    arma::vec wdiag(n);
    arma::vec dinv(n);
    arma::vec ystar(n);
    arma::mat Vstar(n, n);
    arma::mat Vsinv(n, n);
    arma::mat P(n, n);
    arma::mat Minv(m, m, arma::fill::zeros);
    arma::vec sigma_diff(c, arma::fill::zeros);
    arma::vec beta_diff(m, arma::fill::zeros);

    double update_disp = curr_disp;
    double disp_diff = std::numeric_limits<double>::infinity();
    const double disp_tol = 1e-2;
    bool converged = false;
    bool meet_cond = false;
    int iters = 0;

    // broadcast Ginv once per sigma update for the dispersion search
    arma::mat littleG(c, c, arma::fill::zeros);

    while(!meet_cond){
        // ---- working response and weights (both diagonal, kept as vectors) --
        dinv = 1.0 / muvec;
        arma::vec eta = offsets + (X * curr_beta);
        if(qz > 0){
            eta += Z * curr_u.head(qz);
        }
        eta += curr_u.tail(n); // genetic BLUPs enter with an implicit identity design

        ystar = eta + (dinv % (y - muvec));
        wdiag = disp_as_vc ? dinv : ((1.0 / curr_disp) + dinv);

        // The offset is part of the linear predictor but is not a column of X,
        // so it must be removed before the GLS solve. Leaving it in lets the
        // intercept absorb the offset, which then enters eta twice and sends mu
        // to infinity.
        arma::vec ystar_c = ystar - offsets;

        // ---- pseudo-variance and its inverse -------------------------------
        Vstar = buildVstar(wdiag, Z, K, curr_sigma, u_indices);
        if(disp_as_vc){
            Vstar.diag() += sig_a(c);
        }
        bool _vok = arma::inv_sympd(Vsinv, Vstar);
        if(!_vok){
            Rcpp::warning("Pseudovariance is not positive definite - using pseudoinverse");
            Vsinv = arma::pinv(Vstar);
        }

        // ---- REML projection ------------------------------------------------
        arma::mat XtVi = X.t() * Vsinv;
        arma::mat M = XtVi * X;
        bool _mok = arma::inv_sympd(Minv, M);
        if(!_mok){
            Minv = arma::pinv(M);
        }

        if(REML){
            P = Vsinv - (Vsinv * X) * Minv * XtVi;
            P = 0.5 * (P + P.t());
        } else {
            P = Vsinv;
        }

        // ---- variance components -------------------------------------------
        if(!fix_variance){
            arma::vec score_sigma(ctot, arma::fill::zeros);
            arma::mat info_sigma(ctot, ctot, arma::fill::zeros);
            arma::vec Py = P * ystar_c;

            std::vector<arma::mat> PdV(ctot);
            for(int j = 0; j < ctot; j++){
                const arma::mat& dVj = dVa[j];
                PdV[j] = P * dVj;
            }

            for(int j = 0; j < ctot; j++){
                const arma::mat& dVj = dVa[j];
                double lhs = -0.5 * arma::trace(PdV[j]);
                double rhs = 0.5 * arma::as_scalar(Py.t() * dVj * Py);
                score_sigma(j) = lhs + rhs;

                for(int k = j; k < ctot; k++){
                    // trace(A * B) without forming the product
                    double tr = arma::accu(PdV[j] % PdV[k].t());
                    info_sigma(j, k) = 0.5 * tr;
                    if(j != k){
                        info_sigma(k, j) = 0.5 * tr;
                    }
                }
            }

            arma::vec sigma_update = fisherScore(info_sigma, score_sigma, sig_a);

            // The domain of the variance components is [0, Inf). Rather than
            // clamping a negative update straight to the boundary - which pins
            // the component there for every later iteration, because the next
            // Fisher step is taken from the boundary - retreat along the same
            // ascent direction by step-halving until every component is
            // strictly positive. This keeps the search direction and lets a
            // component recover if the data support it.
            arma::vec step = sigma_update - sig_a;
            int halvings = 0;
            while(arma::any((sig_a + step) <= 0.0) && halvings < 30){
                step *= 0.5;
                halvings++;
            }
            sigma_update = sig_a + step;

            // final guard for non-finite or still non-positive components
            for(int j = 0; j < ctot; j++){
                if(!std::isfinite(sigma_update(j)) || sigma_update(j) <= 0.0){
                    sigma_update(j) = constval;
                }
            }

            sigma_diff = arma::abs(sigma_update.head(c) - curr_sigma);
            sig_a = sigma_update;
            curr_sigma = sigma_update.head(c);
            if(disp_as_vc){
                // report the overdispersion on the size scale
                curr_disp = 1.0 / std::max(sig_a(c), 1e-12);
            }
        }

        // ---- fixed effects and BLUPs ---------------------------------------
        // recompute V* with the updated variance components before solving
        Vstar = buildVstar(wdiag, Z, K, curr_sigma, u_indices);
        if(disp_as_vc){
            Vstar.diag() += sig_a(c);
        }
        _vok = arma::inv_sympd(Vsinv, Vstar);
        if(!_vok){
            Vsinv = arma::pinv(Vstar);
        }
        XtVi = X.t() * Vsinv;
        M = XtVi * X;
        _mok = arma::inv_sympd(Minv, M);
        if(!_mok){
            Minv = arma::pinv(M);
        }

        arma::vec beta_update = Minv * (XtVi * ystar_c);
        beta_diff = arma::abs(beta_update - curr_beta);
        curr_beta = beta_update;

        arma::vec resid = ystar_c - (X * curr_beta);
        arma::vec Vresid = Vsinv * resid;

        arma::vec u_update(qz + n, arma::fill::zeros);
        for(int j = 0; j < c - 1; j++){
            arma::uvec jdx = u_indices[j];
            arma::mat Zj = Z.cols(jdx - 1);
            u_update.elem(jdx - 1) = curr_sigma(j) * (Zj.t() * Vresid);
        }
        u_update.tail(n) = curr_sigma(c - 1) * (K * Vresid);
        curr_u = u_update;

        // ---- dispersion ------------------------------------------------------
        if(!fix_variance){
            littleG.zeros();
            for(int j = 0; j < c; j++){
                littleG(j, j) = curr_sigma(j);
            }

            arma::mat Ginv(qz + n, qz + n, arma::fill::zeros);
            for(int j = 0; j < c - 1; j++){
                arma::uvec jdx = u_indices[j];
                for(unsigned int l = 0; l < jdx.n_elem; l++){
                    Ginv(jdx(l) - 1, jdx(l) - 1) = 1.0 / curr_sigma(j);
                }
            }
            Ginv.submat(qz, qz, qz + n - 1, qz + n - 1) = Kinv / curr_sigma(c - 1);

            // The golden-section search resolves the dispersion only to its own
            // tolerance (1e-2), so continuing to update it every iteration
            // injects jitter of that size into W and prevents the fixed effects
            // ever meeting a tighter convergence tolerance. Stop updating once
            // the dispersion has settled, as fitGeneticPLGlmm does.
            if(!fix_dispersion && !disp_as_vc && disp_diff > disp_tol){
                double delta_lo = std::max(1e-2, curr_disp - (curr_disp * 0.5));
                double delta_up = std::min(max_disp, std::max(2e-2, curr_disp * 2.0));
                update_disp = phiGoldenSearch(curr_disp, delta_lo, delta_up, c,
                                              muvec, Ginv, pi, curr_u, curr_sigma, y);
                // phiGoldenSearch maximises the conditional NB likelihood at the
                // current mu, and mu already contains the fitted BLUPs. The
                // random effects have therefore already absorbed the
                // overdispersion, so the search sees almost none left and drives
                // the size parameter towards the Poisson limit. Bounding it stops
                // the runaway; the principled fix is to supply an externally
                // estimated dispersion and set fix_dispersion.
                update_disp = std::min(update_disp, max_disp);
                disp_diff = std::abs(curr_disp - update_disp);
                curr_disp = update_disp;
            }
        }

        // ---- update the linear predictor -------------------------------------
        arma::vec eta_new = offsets + (X * curr_beta);
        if(qz > 0){
            eta_new += Z * curr_u.head(qz);
        }
        eta_new += curr_u.tail(n);
        muvec = arma::exp(eta_new);

        LogicalVector _chk_na = check_na_arma_numeric(muvec);
        if(any(_chk_na).is_true()){
            stop("NA estimates in linear predictor - consider an alternative model");
        }
        LogicalVector _chk_inf = check_inf_arma_numeric(muvec);
        if(any(_chk_inf).is_true()){
            stop("Infinite parameter estimates - consider an alternative model");
        }

        iters++;
        bool _bconv = arma::all(beta_diff < theta_conv);
        bool _sconv = fix_variance ? true : arma::all(sigma_diff < theta_conv);
        bool _ithit = iters >= maxit;

        converged = _bconv && _sconv;
        meet_cond = converged || _ithit;
    }

    // ---- final quantities ----------------------------------------------------
    dinv = 1.0 / muvec;
    arma::vec eta_f = offsets + (X * curr_beta);
    if(qz > 0){
        eta_f += Z * curr_u.head(qz);
    }
    eta_f += curr_u.tail(n);
    ystar = eta_f + (dinv % (y - muvec));
    wdiag = disp_as_vc ? dinv : ((1.0 / curr_disp) + dinv);
    arma::vec ystar_c = ystar - offsets;

    Vstar = buildVstar(wdiag, Z, K, curr_sigma, u_indices);
    if(disp_as_vc){
        Vstar.diag() += sig_a(c);
    }
    bool _vok2 = arma::inv_sympd(Vsinv, Vstar);
    if(!_vok2){
        Vsinv = arma::pinv(Vstar);
    }
    arma::mat XtVi = X.t() * Vsinv;
    arma::mat M = XtVi * X;
    bool _mok2 = arma::inv_sympd(Minv, M);
    if(!_mok2){
        Minv = arma::pinv(M);
    }

    if(REML){
        P = Vsinv - (Vsinv * X) * Minv * XtVi;
        P = 0.5 * (P + P.t());
    } else {
        P = Vsinv;
    }

    arma::vec se = arma::sqrt(Minv.diag());
    arma::vec tscore(m);
    for(int i = 0; i < m; i++){
        tscore(i) = se(i) > 0.0 ? curr_beta(i) / se(i) : NA_REAL;
    }

    littleG.zeros();
    for(int j = 0; j < c; j++){
        littleG(j, j) = curr_sigma(j);
    }
    arma::mat Ginv_f(qz + n, qz + n, arma::fill::zeros);
    for(int j = 0; j < c - 1; j++){
        arma::uvec jdx = u_indices[j];
        for(unsigned int l = 0; l < jdx.n_elem; l++){
            Ginv_f(jdx(l) - 1, jdx(l) - 1) = 1.0 / curr_sigma(j);
        }
    }
    Ginv_f.submat(qz, qz, qz + n - 1, qz + n - 1) = Kinv / curr_sigma(c - 1);
    double loglihood = nbLogLik(muvec, curr_disp, y) - normLogLik(c, Ginv_f, littleG, curr_u, pi);

    // Report W on a common scale for both estimators: the disp_as_vc path
    // carries the constant diagonal as sigma_0 rather than inside wdiag, but
    // 1/curr_disp == sigma_0 by construction, so callers can always rebuild
    // V* as diag(Wdiag) + sum_j sigma_j Zj Zj' + sigma_g K.
    if(disp_as_vc){
        wdiag = dinv + (1.0 / curr_disp);
    }

    arma::mat P_out(1, 1, arma::fill::zeros);
    arma::vec Py_out(1, arma::fill::zeros);
    if(return_projection){
        P_out = P;
        Py_out = P * ystar_c;
    }

    return List::create(_["FE"]=curr_beta, _["RE"]=curr_u, _["Sigma"]=curr_sigma,
                        _["Dispersion"]=curr_disp, _["converged"]=converged, _["Iters"]=iters,
                        _["SE"]=se, _["t"]=tscore, _["P"]=P_out, _["Pystar"]=Py_out,
                        _["ystar"]=ystar, _["Wdiag"]=wdiag, _["VCOV"]=Minv,
                        _["LOGLIHOOD"]=loglihood, _["Vsinv"]=Vsinv);
}


//' Score test for genetic variants against a fitted null model
//'
//' Computes score (Rao) test statistics for a block of genetic variants against
//' a null NB-GLMM fitted by \code{\link{fitGeneticNullGlmm}}. Because the
//' variance components do not depend on the variant being tested, the REML
//' projection \emph{P} and the vector \eqn{P y^*} are computed once per
//' neighbourhood and reused across every variant.
//'
//' For a variant with genotype vector \emph{g} the score statistic is
//' \eqn{U = g' P y^*} with variance \eqn{I = g' P g}, giving
//' \eqn{\chi^2_1 = U^2 / I}. The corresponding one-step effect size estimate is
//' \eqn{\hat{\beta} = U / I} with standard error \eqn{1/\sqrt{I}}, so the Wald
//' statistic formed from these is algebraically identical to the score
//' statistic. Variants are processed as a block so that \eqn{P G} is a single
//' matrix-matrix product rather than a loop of matrix-vector products.
//'
//' @param P mat - the REML projection matrix from the fitted null model.
//' @param Pystar vec - the vector \eqn{P y^*} from the fitted null model.
//' @param G mat - an n x B matrix of genotypes, one column per variant.
//' @param min_variance double - variants whose score variance falls below this
//' value are returned as \code{NA} rather than producing an unstable ratio.
//' This traps monomorphic variants and variants that are collinear with the
//' fixed effects.
//'
//' @details The score test uses the variance components estimated under the
//' null. This is the standard approximation used by genome-wide mixed model
//' methods and is asymptotically equivalent to the Wald test from a full refit;
//' variants passing a screening threshold should be refitted exactly with
//' \code{\link{fitGeneticNullGlmm}} to obtain final effect size estimates.
//'
//' @return A \code{list} with one element per statistic, each a \code{numeric}
//' vector of length B:
//' \describe{
//' \item{\code{Score}:}{the score \eqn{U = g' P y^*}.}
//' \item{\code{Variance}:}{the score variance \eqn{I = g' P g}.}
//' \item{\code{Chisq}:}{the test statistic \eqn{U^2 / I} on 1 degree of freedom.}
//' \item{\code{Beta}:}{the one-step effect size estimate \eqn{U / I}.}
//' \item{\code{SE}:}{the standard error \eqn{1/\sqrt{I}}.}
//' }
//'
//' @author Mike Morgan
//'
//' @examples
//' NULL
//'
//' @name scoreTestGeneticSNPs
//'
// [[Rcpp::export]]
List scoreTestGeneticSNPs(const arma::mat& P, const arma::vec& Pystar,
                          const arma::mat& G, double min_variance = 1e-12){

    const unsigned int n = P.n_rows;
    const unsigned int B = G.n_cols;

    if(G.n_rows != n){
        stop("Genotype matrix rows do not match the dimension of the projection matrix");
    }
    if(Pystar.n_elem != n){
        stop("Pystar length does not match the dimension of the projection matrix");
    }

    // single GEMM for the whole block - this is the only O(n^2 B) term
    arma::mat PG = P * G;

    arma::vec score = G.t() * Pystar;
    arma::vec variance(B);
    arma::vec chisq(B);
    arma::vec beta(B);
    arma::vec se(B);

    for(unsigned int s = 0; s < B; s++){
        double v = arma::dot(G.col(s), PG.col(s));
        if(!std::isfinite(v) || v < min_variance){
            variance(s) = NA_REAL;
            chisq(s) = NA_REAL;
            beta(s) = NA_REAL;
            se(s) = NA_REAL;
            continue;
        }
        variance(s) = v;
        chisq(s) = (score(s) * score(s)) / v;
        beta(s) = score(s) / v;
        se(s) = 1.0 / std::sqrt(v);
    }

    return List::create(_["Score"]=score, _["Variance"]=variance, _["Chisq"]=chisq,
                        _["Beta"]=beta, _["SE"]=se);
}
