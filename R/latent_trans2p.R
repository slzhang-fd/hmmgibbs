
#' @export sample_struc2p
sample_struc2p <- function(C_it_tp, C_it_fp, xcovs, u_tp0, u_fp0, struc_params0){
  # struc_params = list(b_tp, b_fp, d_tp, d_fp, beta_tp, beta_fp, rho_u)
  b_tp <- struc_params0$b_tp
  b_fp <- struc_params0$b_fp
  d_tp <- struc_params0$d_tp
  d_fp <- struc_params0$d_fp
  beta_tp <- struc_params0$beta_tp
  beta_fp <- struc_params0$beta_fp
  rho_u <- struc_params0$rho_u

  N <- length(u_tp0)
  wave_num <- length(C_it_tp) / N
  pp <- nrow(b_tp)

  # sample b_tpfp1
  for(p in 1:pp){
    log_b_tp1_pdf <- function(x){
      coeffs1 <- b_tp[,1]
      coeffs1[p] <- x
      tmp <- xcovs[1:N,] %*% coeffs1 + d_tp * u_tp0
      -1/200 * x^2 + sum(C_it_tp[1:N] * tmp - log(1 + exp(tmp)))
    }
    b_tp[p, 1] <- arms(1, log_b_tp1_pdf, -10, 10, metropolis = FALSE)
    log_b_fp1_pdf <- function(x){
      coeffs1 <- b_fp[,1]
      coeffs1[p] <- x
      tmp <- xcovs[1:N,] %*% coeffs1 + d_fp * u_fp0
      -1/200 * x^2 + sum(C_it_fp[1:N] * tmp - log(1 + exp(tmp)))
    }
    b_fp[p, 1] <- arms(1, log_b_fp1_pdf, -10, 10, metropolis = FALSE)
    # sample b_tpfp_2-t
    for(t in 2:wave_num){
      log_b_tpj_pdf <- function(x){
        coeffst <- b_tp[,t]
        coeffst[p] <- x
        tmp <- xcovs[(t-1)*N+1:N,] %*% coeffst + beta_tp[1] * C_it_tp[1:N + (t-2)*N] + beta_tp[2] * C_it_fp[1:N + (t-2)*N] + u_tp0
        -1/200 * x^2 + sum(C_it_tp[(t-1)*N + (1:N)] * tmp - log(1 + exp(tmp)))
      }
      b_tp[p, t] <- arms(1, log_b_tpj_pdf, -10, 10, metropolis = FALSE)
      log_b_fpj_pdf <- function(x){
        coeffst <- b_fp[,t]
        coeffst[p] <- x
        tmp <- xcovs[(t-1)*N+1:N,] %*% coeffst + beta_fp[1] * C_it_fp[1:N + (t-2)*N] + beta_fp[2] * C_it_tp[1:N + (t-2)*N] + u_fp0
        -1/200 * x^2 + sum(C_it_fp[(t-1)*N + (1:N)] * tmp - log(1 + exp(tmp)))
      }
      b_fp[p, t] <- arms(1, log_b_fpj_pdf, -10, 10, metropolis = FALSE)
    }
  }
  it_inter_tp <- it_inter_fp <- matrix(0, N, wave_num)
  for(t in 1:wave_num){
    it_inter_tp[,t] <- xcovs[(t-1)*N+1:N,] %*% b_tp[,t]
    it_inter_fp[,t] <- xcovs[(t-1)*N+1:N,] %*% b_fp[,t]
  }
  # sample d_tp, d_fp
  log_d_tp_pdf <- function(x){
    tmp <- it_inter_tp[,1] + x * u_tp0
    -1/200 * x^2 + sum(C_it_tp[1:N] * (tmp) - log(1 + exp(tmp)))
  }
  d_tp <- arms(1, log_d_tp_pdf, -10, 10, metropolis = FALSE)
  log_d_fp_pdf <- function(x){
    tmp <- it_inter_fp[,1] + x * u_fp0
    -1/200 * x^2 + sum(C_it_fp[1:N] * (tmp) - log(1 + exp(tmp)))
  }
  d_fp <- arms(1, log_d_fp_pdf, -10, 10, metropolis = FALSE)
  ## sample beta_tp, beta_fp
  log_beta_tp1_pdf <- function(x){
    log_beta_tp1_val <- -1/200 * x^2
    for(t in 2:wave_num){
      tmp <- it_inter_tp[,t] + x * C_it_tp[1:N + (t-2)*N] + beta_tp[2] * C_it_fp[1:N+(t-2)*N] + u_tp0
      log_beta_tp1_val <- log_beta_tp1_val + sum(C_it_tp[(t-1)*N + (1:N)] * tmp - log(1 + exp(tmp)))
    }
    log_beta_tp1_val
  }
  beta_tp[1] <- arms(1, log_beta_tp1_pdf, -10, 10, metropolis = FALSE)
  log_beta_tp2_pdf <- function(x){
    log_beta_tp2_val <- -1/200 * x^2
    for(t in 2:wave_num){
      tmp <- it_inter_tp[,t] + beta_tp[1] * C_it_tp[1:N + (t-2)*N] + x * C_it_fp[1:N+(t-2)*N] + u_tp0
      log_beta_tp2_val <- log_beta_tp2_val + sum(C_it_tp[(t-1)*N + (1:N)] * tmp - log(1 + exp(tmp)))
    }
    log_beta_tp2_val
  }
  beta_tp[2] <- arms(1, log_beta_tp2_pdf, -10, 10, metropolis = FALSE)
  log_beta_fp1_pdf <- function(x){
    log_beta_fp1_val <- -1/200 * x^2
    for(t in 2:wave_num){
      tmp <- it_inter_fp[,t] + x * C_it_fp[1:N + (t-2)*N] + beta_fp[2] * C_it_tp[1:N+(t-2)*N] + u_fp0
      log_beta_fp1_val <- log_beta_fp1_val + sum(C_it_fp[(t-1)*N + (1:N)] * tmp - log(1 + exp(tmp)))
    }
    log_beta_fp1_val
  }
  beta_fp[1] <- arms(1, log_beta_fp1_pdf, -10, 10, metropolis = FALSE)
  log_beta_fp2_pdf <- function(x){
    log_beta_fp2_val <- -1/200 * x^2
    for(t in 2:wave_num){
      tmp <- it_inter_fp[,t] + beta_fp[1] * C_it_fp[1:N + (t-2)*N] + x * C_it_tp[1:N+(t-2)*N] + u_fp0
      log_beta_fp2_val <- log_beta_fp2_val + sum(C_it_fp[(t-1)*N + (1:N)] * tmp - log(1 + exp(tmp)))
    }
    log_beta_fp2_val
  }
  beta_fp[2] <- arms(1, log_beta_fp2_pdf, -10, 10, metropolis = FALSE)
  ## sample rho_u
  log_rho_u_pdf <- function(x){
    -0.5 * N * log(1-x^2) - 0.5 / (1-x^2) * sum(u_tp0^2 - 2*x*u_tp0*u_fp0 + u_fp0^2)
  }
  rho_u <- arms(1, log_rho_u_pdf, -1, 1, metropolis = TRUE, previous = rho_u)
  ## sample u_tp0, u_fp0
  u_tp0 <- sample_u2p_cpp(it_inter_tp, xcovs, d_tp, beta_tp, rho_u, matrix(C_it_tp, N), matrix(C_it_fp, N), u_fp0)
  u_fp0 <- sample_u2p_cpp(it_inter_fp, xcovs, d_fp, beta_fp, rho_u, matrix(C_it_fp, N), matrix(C_it_tp, N), u_tp0)
  # struc_params = list(b_tp, b_fp, d_tp, d_fp, beta_tp, beta_fp, rho_u)
  list(params = list(b_tp = b_tp,
       b_fp = b_fp,
       d_tp = d_tp,
       d_fp = d_fp,
       beta_tp = beta_tp,
       beta_fp = beta_fp,
       rho_u = rho_u),
       u_tp0 = u_tp0,
       u_fp0 = u_fp0,
       it_inter_tp = it_inter_tp,
       it_inter_fp = it_inter_fp)
}
#' @export sample_struc2p
sample_C <- function(y_tpfp, u_tpfp, u_fptp, C_tpfp, C_fptp,
                     alpha_tpfp, alpha_fptp, it_inter_tpfp, it_inter_fptp, d_tpfp,
                     beta_tpfp, beta_fptp){
  N <- length(u_tpfp)
  lambda1 <- t(alpha_tpfp[1,] + outer(alpha_tpfp[2,], rep(1,N)))
  lambda0 <- t(alpha_tpfp[1,] + outer(alpha_tpfp[2,], rep(0,N)))

  tmp <- it_inter_tpfp[,1] + d_tpfp * u_tpfp
  tmp1 <- it_inter_tpfp[,2] + beta_tpfp[1] + beta_tpfp[2] * C_fptp[1:N] + u_tpfp
  tmp2 <- it_inter_fptp[,2] + beta_fptp[1] * C_fptp[1:N] + beta_fptp[2] + u_fptp
  parts1 <- rowSums(lambda1 * y_tpfp[1:N,] - log(1+exp(lambda1))) +
    tmp - log(1+exp(tmp)) +
    C_tpfp[N+1:N] * tmp1 - log(1+exp(tmp1)) +
    C_fptp[N+1:N] * tmp2 - log(1+exp(tmp2))
  tmp10 <- it_inter_tpfp[,2] + beta_tpfp[2] * C_fptp[1:N] + u_tpfp
  tmp20 <- it_inter_fptp[,2] + beta_fptp[1] * C_fptp[1:N] + u_fptp
  parts0 <- rowSums(lambda0 * y_tpfp[1:N,] - log(1+exp(lambda0)))+
    - log(1+exp(tmp)) +
    C_tpfp[N+1:N] * tmp10 - log(1+exp(tmp10)) +
    C_fptp[N+1:N] * tmp20 - log(1+exp(tmp20))
  C_tpfp[1:N] <- rbinom(N, 1, prob = 1 / (1 + exp(parts0 - parts1)))

  tmp <- it_inter_tpfp[,2] + beta_tpfp[1] * C_tpfp[1:N] + beta_tpfp[2] * C_fptp[1:N] + u_tpfp
  tmp1 <- it_inter_tpfp[,3] + beta_tpfp[1] + beta_tpfp[2] * C_fptp[N+1:N] + u_tpfp
  tmp2 <- it_inter_fptp[,3] + beta_fptp[1] * C_fptp[N+1:N] + beta_fptp[2] + u_fptp
  parts1 <- rowSums(lambda1 * y_tpfp[N+1:N,] - log(1+exp(lambda1))) +
    tmp - log(1+exp(tmp)) +
    C_tpfp[2*N+1:N] * tmp1 - log(1+exp(tmp1))+
    C_fptp[2*N+1:N] * tmp2 - log(1+exp(tmp2))
  tmp10 <- it_inter_tpfp[,3] + beta_tpfp[2] * C_fptp[N+1:N] + u_tpfp
  tmp20 <- it_inter_fptp[,3] + beta_fptp[1] * C_fptp[N+1:N] + u_fptp
  parts0 <- rowSums(lambda0 * y_tpfp[N+1:N,] - log(1+exp(lambda0))) +
    - log(1+exp(tmp)) +
    C_tpfp[2*N+1:N] * tmp10 - log(1+exp(tmp10))+
    C_fptp[2*N+1:N] * tmp20 - log(1+exp(tmp20))
  C_tpfp[N+1:N] <- rbinom(N, 1, prob = 1 / (1 + exp(parts0 - parts1)))

  tmp <- it_inter_tpfp[,3] + beta_tpfp[1] * C_tpfp[N+1:N] + beta_tpfp[2] * C_fptp[N+1:N] + u_tpfp
  parts1 <- rowSums(lambda1 * y_tpfp[2*N+1:N,] - log(1+exp(lambda1))) + tmp - log(1+exp(tmp))
  parts0 <- rowSums(lambda0 * y_tpfp[2*N+1:N,] - log(1+exp(lambda0))) - log(1+exp(tmp))
  C_tpfp[2*N+1:N] <- rbinom(N, 1, prob = 1 / (1 + exp(parts0 - parts1)))
  C_tpfp
}

#' @export latent_trans2p
latent_trans2p <- function(y_tp, y_fp, xcovs,
                           C_it_tp, C_it_fp, u_tp0, u_fp0,
                           alpha_tp0, alpha_fp0, struc_params0,
                           mcmc_len){
  struc_params_draw <- matrix(0, mcmc_len, length(unlist(struc_params0)))
  alpha_tp_draws <- matrix(0, mcmc_len, 2*J)
  alpha_fp_draws <- matrix(0, mcmc_len, 2*J)

  for(iter in 1:mcmc_len){
    cat('iter= ', iter, '\n')
    ## sample measurement parameters: alpha_tp, alpha_fp
    alpha_tp0 <- sample_alpha_cpp(y_tp, C_it_tp, alpha_tp0)
    alpha_fp0 <- sample_alpha_cpp(y_fp, C_it_fp, alpha_fp0)

    # sample structure parameters:
    struc_temp <- sample_struc2p(C_it_tp, C_it_fp, xcovs, u_tp0, u_fp0, struc_params0)
    struc_params0 <- struc_temp$params
    u_tp0 <- struc_temp$u_tp0
    u_fp0 <- struc_temp$u_fp0
    ## sample C_it_tp
    C_it_tp <- sample_C(y_tp, u_tp0, u_fp0, C_it_tp, C_it_fp, alpha_tp0, alpha_fp0,
                        struc_params0$it_inter_tp, struc_params0$it_inter_fp, struc_params0$d_tp,
                        struc_params0$beta_tp, struc_params0$beta_fp)
    C_it_fp <- sample_C(y_fp, u_fp0, u_tp0, C_it_fp, C_it_tp, alpha_fp0, alpha_tp0,
                        struc_params0$it_inter_fp, struc_params0$it_inter_tp, struc_params0$d_fp,
                        struc_params0$beta_fp, struc_params0$beta_tp)
    # store results
    alpha_tp_draws[iter,] <- c(alpha_tp0)
    alpha_fp_draws[iter,] <- c(alpha_fp0)
    struc_params_draw[iter,] <- unlist(struc_params0)
  }
  list(alpha_tp_draws=alpha_tp_draws,
       alpha_fp_draws=alpha_fp_draws,
       struc_params_draws=struc_params_draw)
}
