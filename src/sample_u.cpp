// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "arms_ori.h"

struct log_u_param{
  arma::vec bj;
  double d1;
  double beta;
  double sigma2;
  arma::vec Ci_tp;
};
double log_u_pdf(double x, void* params){
  struct log_u_param *d;
  d = static_cast<struct log_u_param *> (params);

  int wave_num = d->bj.n_rows;
  double tmp = d->bj(0) + d->d1 * x;
  double log_u_val = -0.5 / d->sigma2 * x * x + d->Ci_tp(0) * tmp - std::log(1.0 + std::exp(tmp));
  for(unsigned int t=1; t<wave_num; ++t){
    tmp = d->bj(t) + d->beta * d->Ci_tp(t-1) + x;
    log_u_val += d->Ci_tp(t) * tmp - std::log(1.0 + std::exp(tmp));
  }
  return log_u_val;
}
//' @export
// [[Rcpp::export]]
double log_u_pdf_check(double x, arma::vec bj, double d1, double beta, double sigma2, arma::vec Ci_tp){
  int wave_num = bj.n_rows;
  double tmp = bj(0) + d1 * x;
  double log_u_val = -0.5 / sigma2 * x * x + Ci_tp(0) * tmp - std::log(1.0 + std::exp(tmp));
  for(unsigned int t=1; t<wave_num; ++t){
    tmp = bj(t) + beta * Ci_tp(t-1) + x;
    log_u_val += Ci_tp(t) * tmp - std::log(1.0 + std::exp(tmp));
  }
  return log_u_val;
}
//' @export
// [[Rcpp::export]]
arma::vec sample_u_cpp_check(arma::vec bj, double d1, double beta, double sigma2,
                       arma::mat C_mat, int nsamp, int ii){
  int err, ninit = 4, npoint = 100, ncent = 0;
  int neval;
  double xl = -100.0, xr = 100.0;
  double xinit[10]={-10.0, -2.0, 2.0, 10.0}, xsamp[nsamp], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  double xprev = 0.0;

  arma::vec res(nsamp);

  log_u_param log_u_data;
  log_u_data.bj = bj;
  log_u_data.d1 = d1;
  log_u_data.beta = beta;
  log_u_data.sigma2 = sigma2;

  log_u_data.Ci_tp = C_mat.row(ii).t();
  err = arms(xinit,ninit,&xl,&xr,log_u_pdf,&log_u_data,&convex,
             npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
  if(err>0){
    Rprintf("error code: %d", err);
    Rcpp::stop("\n");
  }
  for(int i=0;i<nsamp;++i){
    res(i) = xsamp[i];
  }

  return res;
}
// [[Rcpp::export]]
arma::vec sample_u_cpp(arma::vec bj, double d1, double beta, double sigma2,
                       arma::mat C_mat){
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 0;
  int neval;
  double xl = -100.0, xr = 100.0;
  double xinit[10]={-10.0, -2.0, 2.0, 10.0}, xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  double xprev = 0.0;

  unsigned int N = C_mat.n_rows;
  arma::vec u_tp_res(N);

  log_u_param log_u_data;
  log_u_data.bj = bj;
  log_u_data.d1 = d1;
  log_u_data.beta = beta;
  log_u_data.sigma2 = sigma2;

  for(unsigned int i=0;i<N;++i){
    log_u_data.Ci_tp = C_mat.row(i).t();
    err = arms(xinit,ninit,&xl,&xr,log_u_pdf,&log_u_data,&convex,
               npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
    if(err>0){
      Rprintf("error code: %d", err);
      Rcpp::stop("\n");
    }
    u_tp_res(i) = xsamp[0];
  }
  return u_tp_res;
}
