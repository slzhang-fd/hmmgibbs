// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "arms_ori.h"

// log_alpha0_j <- function(x){
//   tmp <- x + alpha_tp0[2,j] * C_it_tp
//   -1/200 * x^2 + sum(y_tp[,j] * tmp - log(1 + exp(tmp)))
// }
// log_alpha1_j <- function(x){
//   tmp <- alpha_tp0[1,j] + x * C_it_tp
//   -1/200 * x^2 + sum(y_tp[,j] * tmp - log(1 + exp(tmp)))
// }
struct log_alpha_param{
  arma::vec y_tp_fp_j;
  arma::vec C_tp_fp;
  double alpha10_j;
};
double log_alpha0_pdf(double x, void* params){
  struct log_alpha_param *d;
  d = static_cast<struct log_alpha_param *> (params);

  arma::vec tmp = x + d->alpha10_j * d->C_tp_fp;
  return -1.0/200 * x * x +
    arma::accu(d->y_tp_fp_j % tmp - arma::log(1.0 + arma::exp(tmp)));
}
double log_alpha1_pdf(double x, void* params){
  struct log_alpha_param *d;
  d = static_cast<struct log_alpha_param *> (params);

  arma::vec tmp = d->alpha10_j + x * d->C_tp_fp;
  return -1.0/200 * x * x +
    arma::accu(d->y_tp_fp_j % tmp - arma::log(1.0 + arma::exp(tmp)));
}
//' @export
// [[Rcpp::export]]
arma::mat sample_alpha_cpp(arma::mat y_tpfp, arma::vec C_tpfp,
                           arma::mat alpha_tpfp){
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 0;
  int neval;
  double xl = -100.0, xr = 100.0;
  double xinit[10]={-10.0, -2.0, 2.0, 10.0}, xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  double xprev = 0.0;

  unsigned int N = y_tpfp.n_rows;
  unsigned int J = y_tpfp.n_cols;

  log_alpha_param log_alpha_data;
  log_alpha_data.C_tp_fp = C_tpfp;
  for(unsigned int j = 0; j<J;++j){
    // sample alpha_tp_0j
    // Rcpp::Rcout << "j= " << j << std::endl;
    log_alpha_data.y_tp_fp_j = y_tpfp.col(j);
    log_alpha_data.alpha10_j = alpha_tpfp(1, j);
    err = arms(xinit,ninit,&xl,&xr,log_alpha0_pdf,&log_alpha_data,&convex,
               npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
    if(err>0){
      Rprintf("error code: %d", err);
      Rcpp::stop("\n");
    }
    alpha_tpfp(0, j) = xsamp[0];
    // sample alpha_tp_1j
    log_alpha_data.alpha10_j = alpha_tpfp(0, j);
    err = arms(xinit,ninit,&xl,&xr,log_alpha1_pdf,&log_alpha_data,&convex,
               npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
    if(err>0){
      Rprintf("error code: %d", err);
      Rcpp::stop("\n");
    }
    alpha_tpfp(1, j) = xsamp[0];
  }
  return alpha_tpfp;
}
