// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sample_alpha_cpp
arma::mat sample_alpha_cpp(arma::mat y_tpfp, arma::vec C_tpfp, arma::mat alpha_tpfp);
RcppExport SEXP _hmmgibbs_sample_alpha_cpp(SEXP y_tpfpSEXP, SEXP C_tpfpSEXP, SEXP alpha_tpfpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y_tpfp(y_tpfpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type C_tpfp(C_tpfpSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_tpfp(alpha_tpfpSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_alpha_cpp(y_tpfp, C_tpfp, alpha_tpfp));
    return rcpp_result_gen;
END_RCPP
}
// log_u_pdf_check
double log_u_pdf_check(double x, arma::vec bj, double d1, double beta, double sigma2, arma::vec Ci_tp);
RcppExport SEXP _hmmgibbs_log_u_pdf_check(SEXP xSEXP, SEXP bjSEXP, SEXP d1SEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP Ci_tpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bj(bjSEXP);
    Rcpp::traits::input_parameter< double >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ci_tp(Ci_tpSEXP);
    rcpp_result_gen = Rcpp::wrap(log_u_pdf_check(x, bj, d1, beta, sigma2, Ci_tp));
    return rcpp_result_gen;
END_RCPP
}
// sample_u_cpp_check
arma::vec sample_u_cpp_check(arma::vec bj, double d1, double beta, double sigma2, arma::mat C_mat, int nsamp, int ii);
RcppExport SEXP _hmmgibbs_sample_u_cpp_check(SEXP bjSEXP, SEXP d1SEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP C_matSEXP, SEXP nsampSEXP, SEXP iiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type bj(bjSEXP);
    Rcpp::traits::input_parameter< double >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_mat(C_matSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_u_cpp_check(bj, d1, beta, sigma2, C_mat, nsamp, ii));
    return rcpp_result_gen;
END_RCPP
}
// sample_u_cpp
arma::vec sample_u_cpp(arma::vec bj, double d1, double beta, double sigma2, arma::mat C_mat);
RcppExport SEXP _hmmgibbs_sample_u_cpp(SEXP bjSEXP, SEXP d1SEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP C_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type bj(bjSEXP);
    Rcpp::traits::input_parameter< double >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_mat(C_matSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_u_cpp(bj, d1, beta, sigma2, C_mat));
    return rcpp_result_gen;
END_RCPP
}
// sample_u2p_cpp
arma::vec sample_u2p_cpp(arma::mat it_inter, arma::mat xcovs, double d_tpfp, arma::vec beta_tpfp, double rho_u, arma::mat C_tpfp_mat, arma::mat C_fptp_mat, arma::vec u_fp_tp);
RcppExport SEXP _hmmgibbs_sample_u2p_cpp(SEXP it_interSEXP, SEXP xcovsSEXP, SEXP d_tpfpSEXP, SEXP beta_tpfpSEXP, SEXP rho_uSEXP, SEXP C_tpfp_matSEXP, SEXP C_fptp_matSEXP, SEXP u_fp_tpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type it_inter(it_interSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xcovs(xcovsSEXP);
    Rcpp::traits::input_parameter< double >::type d_tpfp(d_tpfpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_tpfp(beta_tpfpSEXP);
    Rcpp::traits::input_parameter< double >::type rho_u(rho_uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_tpfp_mat(C_tpfp_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_fptp_mat(C_fptp_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u_fp_tp(u_fp_tpSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_u2p_cpp(it_inter, xcovs, d_tpfp, beta_tpfp, rho_u, C_tpfp_mat, C_fptp_mat, u_fp_tp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hmmgibbs_sample_alpha_cpp", (DL_FUNC) &_hmmgibbs_sample_alpha_cpp, 3},
    {"_hmmgibbs_log_u_pdf_check", (DL_FUNC) &_hmmgibbs_log_u_pdf_check, 6},
    {"_hmmgibbs_sample_u_cpp_check", (DL_FUNC) &_hmmgibbs_sample_u_cpp_check, 7},
    {"_hmmgibbs_sample_u_cpp", (DL_FUNC) &_hmmgibbs_sample_u_cpp, 5},
    {"_hmmgibbs_sample_u2p_cpp", (DL_FUNC) &_hmmgibbs_sample_u2p_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_hmmgibbs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
