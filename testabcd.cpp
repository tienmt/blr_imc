// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::export]]
SEXP mymatmul(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> D){
  Eigen::MatrixXd C = A * B * D;
  return Rcpp::wrap(C);
}