#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Compute pseudovariance partial derivatives
//'
//' Compute the partial derivatives of the pseudovariance as t(Z[, i]) %*% Z[, i]
//' for the ith variance component
//'
//' @param Matrix x - the fully expanded Z matrix that maps observations to
//' random effect variables
//' @param List rlevels - a list that maps the random effect variable to the
//' individual levels
//' @param List dimnames - a list of the matrix `x` dimension names.
// [[Rcpp::export]]
List pseudovarPartial(arma::mat x, List rlevels, StringVector cnames){
    // this currently doesn't support sparse matrices - it's not super clear how to do
    // that concretely without defining some sparse matrix class somewhere along the line

    int items = rlevels.size();
    List outlist(items);

    for(int i = 0; i < items; i++){
        StringVector lelements = rlevels[i];
        IntegerVector icol = match(lelements, cnames); // need  to add a check here in case of NAs or no matches
        int low = min(icol)-1; // need indexing correction from R to C++
        int hi = max(icol)-1;

        // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
        arma::sp_mat omat(x.n_rows, x.n_rows);
        omat = x.cols(low, hi) * x.cols(low, hi).t();
        outlist[i] = omat;
    }

    return outlist;
}
