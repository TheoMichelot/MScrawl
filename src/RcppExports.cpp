// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// kalman_rcpp
double kalman_rcpp(arma::mat& data, arma::vec param, arma::mat& Hmat, arma::vec a0, arma::mat P0);
RcppExport SEXP _MScrawl_kalman_rcpp(SEXP dataSEXP, SEXP paramSEXP, SEXP HmatSEXP, SEXP a0SEXP, SEXP P0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    rcpp_result_gen = Rcpp::wrap(kalman_rcpp(data, param, Hmat, a0, P0));
    return rcpp_result_gen;
END_RCPP
}
// smooth_rcpp
Rcpp::List smooth_rcpp(arma::mat& data, arma::vec param, arma::mat& Hmat, arma::vec a0, arma::mat P0);
RcppExport SEXP _MScrawl_smooth_rcpp(SEXP dataSEXP, SEXP paramSEXP, SEXP HmatSEXP, SEXP a0SEXP, SEXP P0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_rcpp(data, param, Hmat, a0, P0));
    return rcpp_result_gen;
END_RCPP
}
