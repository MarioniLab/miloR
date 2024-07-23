#include <RcppArmadillo.h>
#include <RcppMLCommon.hpp>
#include <RcppEigen.h>
#include <RcppML/nmf.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppML)]]

// solve HE as a constrained QP using RcppML::nnls
arma::vec solveQP(arma::mat A_arma, arma::vec b_arma, arma::vec x_arma){
    // solve the constrained QP problem Ax = y
    // s.t. x >= 0
    Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(A_arma.memptr(), A_arma.n_rows, A_arma.n_cols);
    Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_arma.memptr(), b_arma.n_elem, 1); // column vector
    Eigen::MatrixXd x = Eigen::Map<Eigen::VectorXd>(x_arma.memptr(), x_arma.n_elem, 1); // column vector

    // c_nnls works in-place on an Eigen matrix
    c_nnls(A, b, x, 0);

    // convert x back to an arma::vec
    arma::vec x_res(x.data(), x.size());
    return x_res;
}
