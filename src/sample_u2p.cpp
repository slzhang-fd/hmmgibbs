// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "arms_ori.h"

struct log_u2p_param{
  arma::vec b_tpfp;
  double d_tpfp;
  arma::vec beta_tp_fp;
  arma::vec Ci_tpfp;
  arma::vec Ci_fptp;
  double ui_fptp;
  double rho_u;
};
double log_u2p_pdf(double x, void* params){
  struct log_u2p_param *d;
  d = static_cast<struct log_u2p_param *> (params);

  int wave_num = d->b_tpfp.n_rows;
  double tmp = d->b_tpfp(0) + d->d_tpfp * x;
  double log_u_val = -0.5 / (1-d->rho_u*d->rho_u) * (x * x - 2 * d->rho_u * x * d->ui_fptp);
  log_u_val += d->Ci_tpfp(0) * tmp - std::log(1.0 + std::exp(tmp));
  for(unsigned int t=1; t<wave_num; ++t){
    tmp = d->b_tpfp(t) + d->beta_tp_fp(0) * d->Ci_tpfp(t-1) + d->beta_tp_fp(1) * d->Ci_fptp(t-1) + x;
    log_u_val += d->Ci_tpfp(t) * tmp - std::log(1.0 + std::exp(tmp));
  }
  return log_u_val;
}

// [[Rcpp::export]]
arma::vec sample_u2p_cpp(arma::vec b_tpfp, double d_tpfp,
                         arma::vec beta_tpfp, double rho_u,
                         arma::mat C_tpfp_mat, arma::mat C_fptp_mat,
                         arma::vec u_fp_tp){
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 0;
  int neval;
  double xl = -100.0, xr = 100.0;
  double xinit[10]={-10.0, -2.0, 2.0, 10.0}, xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 1;
  double xprev = 0.0;

  unsigned int N = C_tpfp_mat.n_rows;
  arma::vec u_tpfp_res(N);

  log_u2p_param log_u2p_data;
  log_u2p_data.rho_u = rho_u;
  log_u2p_data.b_tpfp = b_tpfp;
  log_u2p_data.beta_tp_fp = beta_tpfp;
  log_u2p_data.d_tpfp = d_tpfp;
  for(unsigned int i=0;i<N;++i){
    log_u2p_data.Ci_tpfp = C_tpfp_mat.row(i).t();
    log_u2p_data.Ci_fptp = C_fptp_mat.row(i).t();
    log_u2p_data.ui_fptp = u_fp_tp(i);
    err = arms(xinit,ninit,&xl,&xr,log_u2p_pdf,&log_u2p_data,&convex,
               npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
    if(err>0){
      Rprintf("error code: %d", err);
      Rcpp::stop("\n");
    }
    u_tpfp_res(i,0) = xsamp[0];
  }
  return u_tpfp_res;
}
