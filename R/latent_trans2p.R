#' @export latent_trans2p
latent_trans2p <- function(y_tp, y_fp, C_it_tp, C_it_fp,
                           u_tp0, u_fp0, alpha_tp0, alpha_fp0,
                           struc_params0, free_params, mcmc_len){
  bj <- struc_params0[1:wave_num]
  d1 <- struc_params0[wave_num+1]
  beta <- struc_params0[wave_num+2]
  sigma2 <- struc_params0[wave_num+3]

  struc_params_draw <- matrix(0, mcmc_len, length(struc_params0))
  alpha_tp_draws <- matrix(0, mcmc_len, 2*J)
  alpha_fp_draws <- matrix(0, mcmc_len, 2*J)

  for(iter in 1:mcmc_len){
    cat('iter= ', iter, '\n')
    for(j in 1:J){
      log_alpha0_j <- function(x){
        tmp <- x + alpha_tp0[2,j] * C_it_tp
        -1/200 * x^2 + sum(y_tp[,j] * tmp - log(1 + exp(tmp)))
      }
      log_alpha1_j <- function(x){
        tmp <- alpha_tp0[1,j] + x * C_it_tp
        -1/200 * x^2 + sum(y_tp[,j] * tmp - log(1 + exp(tmp)))
      }
      alpha_tp0[1,j] <- arms(1, log_alpha0_j, -10, 10, metropolis = FALSE)
      alpha_tp0[2,j] <- arms(1, log_alpha1_j, -10, 10, metropolis = FALSE)
    }
    # sample b1
    if(free_params[1]){
      log_b1_pdf <- function(x){
        tmp <- x + d1 * u_tp0
        -1/200 * x^2 + sum(C_it_tp[1:N] * (tmp) - log(1 + exp(tmp)))
      }
      bj[1] <- arms(1, log_b1_pdf, -10, 10, metropolis = FALSE)
    }
    # sample b2-bt
    for(t in 2:wave_num){
      log_bj_pdf <- function(x){
        tmp <- x + beta * C_it_tp[1:N + (t-2)*N] + u_tp0
        -1/200 * x^2 + sum(C_it_tp[(t-1)*N + (1:N)] * (tmp) - log(1 + exp(tmp)))
      }
      if(free_params[t]){
        bj[t] <- arms(1, log_bj_pdf, -10, 10, metropolis = FALSE)
      }
    }
    # sample d1
    if(free_params[4]){
      log_d1_pdf <- function(x){
        tmp <- bj[1] + x * u_tp0
        -1/200 * x^2 + sum(C_it_tp[1:N] * (tmp) - log(1 + exp(tmp)))
      }
      d1 <- arms(1, log_d1_pdf, -10, 10, metropolis = FALSE)
    }
    ## sample beta
    if(free_params[5]){
      log_beta_pdf <- function(x){
        log_beta_val <- -1/200 * x^2
        for(t in 2:wave_num){
          tmp <- bj[t] + x * C_it_tp[1:N + (t-2)*N] + u_tp0
          log_beta_val <- log_beta_val + sum(C_it_tp[(t-1)*N + (1:N)] * (tmp) - log(1 + exp(tmp)))
        }
        log_beta_val
      }
      beta <- arms(1, log_beta_pdf, -10, 10, metropolis = FALSE)
    }
    # ## sample sigma2
    if(free_params[6]){
      sigma2 <- rinvgamma(1, 2 + N/2, 1 + sum(u_tp0^2)/2)
    }
    # ## sample u_tp0
    u_tp0 <- sample_u_cpp(bj, d1, beta, sigma2, matrix(C_it_tp, N))
    # for(i in 1:N){
    #   log_u_pdf <- function(x){
    #     tmp <- bj[1] + d1 * x
    #     log_u_val <- -0.5 / sigma2 * x^2 + C_it_tp[i] * tmp - log(1+exp(tmp))
    #     for(t in 2:wave_num){
    #       tmp <- bj[t] + beta * C_it_tp[(t-2)*N + i] + x
    #       log_u_val <- log_u_val + C_it_tp[(t-1)*N+i] * tmp - log(1+exp(tmp))
    #     }
    #     log_u_val
    #   }
    #   u_tp0[i] <- arms(1, log_u_pdf, -10, 10, metropolis = FALSE)
    # }
    ## sample C_it_tp
    lambda1 <- t(alpha_tp0[1,] + outer(alpha_tp0[2,], rep(1,N)))
    lambda0 <- t(alpha_tp0[1,] + outer(alpha_tp0[2,], rep(0,N)))

    parts1 <- rowSums(lambda1 * y_tp[1:N,] - log(1+exp(lambda1)))
    tmp <- bj[1] + d1 * u_tp0
    tmp1 <- bj[2] + beta + u_tp0
    parts1 <- parts1 + tmp - log(1+exp(tmp)) + C_it_tp[N+1:N] * tmp1 - log(1+exp(tmp1))
    parts0 <- rowSums(lambda0 * y_tp[1:N,] - log(1+exp(lambda0)))
    tmp0 <- bj[2] + u_tp0
    parts0 <- parts0 - log(1+exp(tmp)) + C_it_tp[N+1:N] * tmp0 - log(1+exp(tmp0))
    C_it_tp[1:N] <- rbinom(N, 1, prob = 1 / (1 + exp(parts0 - parts1)))

    parts1 <- rowSums(lambda1 * y_tp[N+1:N,] - log(1+exp(lambda1)))
    tmp <- bj[2] + beta * C_it_tp[1:N] + u_tp0
    tmp1 <- bj[3] + beta + u_tp0
    parts1 <- parts1 + tmp - log(1+exp(tmp)) + C_it_tp[2*N+1:N] * tmp1 - log(1+exp(tmp1))
    parts0 <- rowSums(lambda0 * y_tp[N+1:N,] - log(1+exp(lambda0)))
    tmp0 <- bj[3] + u_tp0
    parts0 <- parts0 - log(1+exp(tmp)) + C_it_tp[2*N+1:N] * tmp0 - log(1+exp(tmp0))
    C_it_tp[N+1:N] <- rbinom(N, 1, prob = 1 / (1 + exp(parts0 - parts1)))

    parts1 <- rowSums(lambda1 * y_tp[2*N+1:N,] - log(1+exp(lambda1)))
    tmp <- bj[3] + beta * C_it_tp[N+1:N] + u_tp0
    parts1 <- parts1 + tmp - log(1+exp(tmp))
    parts0 <- rowSums(lambda0 * y_tp[2*N+1:N,] - log(1+exp(lambda0)))
    parts0 <- parts0 - log(1+exp(tmp))
    C_it_tp[2*N+1:N] <- rbinom(N, 1, prob = 1 / (1 + exp(parts0 - parts1)))
    #
    struc_params_draw[iter,] <- c(bj, d1, beta, sigma2)
    alpha_tp_draws[iter,] <- c(alpha_tp0)
  }
  list(alpha_tp_draws, struc_params_draw)
}
