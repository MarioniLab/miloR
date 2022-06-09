#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "pseudovarPartial.h"
using namespace Rcpp;

//' Compute pseudovariance partial derivatives
//'
//' Compute the partial derivatives of the pseudovariance as t(Z[, i]) %*% Z[, i]
//' for the ith variance component
//'
//' @param matrix x - the fully expanded Z matrix that maps observations to
//' random effect variables
//' @param list rlevels - a list that maps the random effect variable to the
//' individual levels
//' @param vec cnames - a StringVector of column names from fully expanded Z matrix
// [[Rcpp::export]]
List pseudovarPartial(arma::mat x, List rlevels, StringVector cnames){
    // this currently doesn't support sparse matrices - it's not super clear how to do
    // that concretely without defining some sparse matrix class somewhere along the line

    unsigned int items = rlevels.size();
    List outlist(items);

    for(unsigned int i = 0; i < items; i++){
        StringVector lelements = rlevels[i];
        IntegerVector icol = match(lelements, cnames); // need  to add a check here in case of NAs or no matches
        int low = min(icol)-1; // need indexing correction from R to C++
        int hi = max(icol)-1;

        // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
        arma::sp_mat omat(x.cols(low, hi) * x.cols(low, hi).t());
        outlist[i] = omat;
    }

    return outlist;
}


// [[Rcpp::export]]
List pseudovarPartial_C(arma::mat Z, List u_indices){
    // A Rcpp specific implementation that uses positional indexing rather than character indexes
    unsigned int items = u_indices.size();
    List outlist(items);

    for(unsigned int i = 0; i < items; i++){
        arma::uvec icols = u_indices[i];

        // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
        arma::mat omat(Z.cols(icols - 1) * Z.cols(icols - 1).t());
        outlist[i] = omat;
    }

    return outlist;
}


List pseudovarPartial_P(List V_partial, const arma::mat& P){
    // A Rcpp specific implementation that uses positional indexing rather than character indexes
    unsigned int items = V_partial.size();
    List outlist(items);

    for(unsigned int i = 0; i < items; i++){
        // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
        arma::mat _omat = V_partial(i);
        arma::mat omat(P * _omat);
        outlist[i] = omat;
    }

    return outlist;

}


List pseudovarPartial_G(arma::mat Z, const arma::mat& K, List u_indices){
    // A Rcpp specific implementation that uses positional indexing rather than character indexes
    unsigned int items = u_indices.size();
    List outlist(items);

    for(unsigned int i = 0; i < items; i++){
        if(i == items - 1){
            outlist(i) = K; // K is equivalent to ZZ^T
        } else{
            arma::uvec icols = u_indices[i];

            // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
            arma::mat _omat(Z.cols(icols - 1) * Z.cols(icols - 1).t());
            outlist[i] = _omat;
        }
    }

    return outlist;

}
