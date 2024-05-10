#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "pseudovarPartial.h"
using namespace Rcpp;

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
        arma::mat omat(x.cols(low, hi) * x.cols(low, hi).t());
        outlist[i] = omat;
    }

    return outlist;
}


List pseudovarPartial_C(arma::mat Z, List u_indices){
    // A Rcpp specific implementation that uses positional indexing rather than character indexes
    unsigned int items = u_indices.size();
    List outlist(items);

    for(unsigned int i = 0; i < items; i++){
        arma::uvec icols = u_indices[i];

        arma::mat omat(Z.cols(icols - 1) * Z.cols(icols - 1).t());
        outlist[i] = omat;
    }

    return outlist;
}


List pseudovarPartial_P(List V_partial, const arma::mat& P){
    // A Rcpp specific implementation that uses positional indexing rather than character indexes
    // don't be tempted to sparsify this - the overhead of casting is too expensive
    unsigned int items = V_partial.size();
    unsigned int n = P.n_cols;
    List outlist(items);
    double temp_value;

    for(unsigned int i = 0; i < items; i++){
        // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
        arma::mat _omat = V_partial(i);
        arma::mat omat(n, n);

        // Can we turn this into a for loop and use OpenMP?
        for (int j = 0; j < n; j++) {
            // j = rows of P
            #pragma omp parallel for reduction(+:temp_value)
            for (int k = 0; k < n; k++) {
                // k = columns of P
                temp_value = 0.0;
                for(int a=0; a < n; a++){
                    temp_value += P(j, a) * _omat(a, k);
                }
                omat(j, k) = temp_value; // Apply P again for symmetry
            }
        }

        outlist[i] = omat;
    }

    return outlist;

}


List pseudovarPartial_V(List V_partial, const arma::mat& V_star_inv){
    // A Rcpp specific implementation that uses positional indexing rather than character indexes
    // don't be tempted to sparsify this - the overhead of casting is too expensive
    unsigned int items = V_partial.size();
    List outlist(items);

    for(unsigned int i = 0; i < items; i++){
        // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
        arma::mat _omat = V_partial(i);
        arma::mat omat(V_star_inv * _omat);
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
            arma::uvec icols = u_indices[i];
            arma::mat _omat(Z.cols(icols - 1) * K * Z.cols(icols - 1).t());
            outlist[i] = _omat; // K is equivalent to ZZ^T
        } else{
            arma::uvec icols = u_indices[i];
            // Need to output an S4 object - arma::sp_mat uses implicit interconversion for support dg Matrices
            arma::mat _omat(Z.cols(icols - 1) * Z.cols(icols - 1).t());
            outlist[i] = _omat;
        }
    }

    return outlist;

}
