library(armspp)
library(MCMCpack)
library(hmmgibbs)
library(MASS)

##################################
#### Simulation parameters setting & data generation
### No covariates in the current simulation
### No missing responses either
##################################

N <- 16439
wave_num <- 3
J <- 8

## coefficients for pi
alpha_tp <- rbind(rnorm(J, 0,1), runif(J))
alpha_fp <- rbind(rnorm(J, 0,1), runif(J))

## random effect u_tp, u_fp
sig2_u_tp <- 1
sig2_u_fp <- 1
rho_u <- 0.5
Sigma_u <- matrix(0, 2,2)
Sigma_u[1,1] <- sig2_u_tp
Sigma_u[2,2] <- sig2_u_fp
Sigma_u[1,2] <- Sigma_u[2,1] <- rho_u * sqrt(sig2_u_tp * sig2_u_fp)
u_tp_fp <- mvrnorm(N, rep(0, 2), Sigma_u)

## time varying covariates
pp <- 3
xcovs <- cbind(rep(1,N*wave_num),matrix(rnorm(N*wave_num*(pp-1)), nrow=N*wave_num))

## time variant intercepts and coefficients for C (random groups)
# b_tp <- c(0, 1, 2)
# b_fp <- c(1, 0, -1)
b_tp <- matrix(rnorm(pp*wave_num, sd = 0.5),nrow = pp)
b_fp <- matrix(rnorm(pp*wave_num, sd = 0.5),nrow = pp)

## coefficients for cross-lagged effect
beta_tp <- c(3, -1)
beta_fp <- c(2, -0.8)

## scale parameter for random effects in wave 1
d_tp <- 1
d_fp <- 1.5


y_tp <- y_fp <- matrix(0, N*wave_num, J)
tind <- rep(1:3,N)
ind <- rep(1:N, each=3)

## generate responses for wave 1
## generate lambda_1^{tp}
tmp <- xcovs %*% b_tp[,1] + d_tp * u_tp_fp[,1]
lambda_1_tp <- 1 / (1+exp(-tmp))
tmp <- xcovs %*% b_fp[,1] + d_fp * u_tp_fp[,2]
lambda_1_fp <- 1 / (1+exp(-tmp))

## generate latent class variables C_1_tp, C_1_fp
C_1_tp <- rbinom(N, 1, prob = lambda_1_tp)
C_1_fp <- rbinom(N, 1, prob = lambda_1_fp)

for(i in 1:N){
  pi_tp_tmp <- alpha_tp[1,] + C_1_tp[i] * alpha_tp[2,]
  y_tp[i,] <- rbinom(J, 1, prob = 1 / (1 + exp(-pi_tp_tmp)))
  pi_fp_tmp <- alpha_fp[1,] + C_1_fp[i] * alpha_fp[2,]
  y_fp[i,] <- rbinom(J, 1, prob = 1 / (1 + exp(-pi_fp_tmp)))
}
## generate responses for wave 2
tmp <- xcovs %*% b_tp[,2] + beta_tp[1] * C_1_tp + beta_tp[2] * C_1_fp + u_tp_fp[,1]
lambda_2_tp <- 1 / (1 + exp(-tmp))
tmp <- xcovs %*% b_fp[,2] + beta_fp[1] * C_1_fp + beta_fp[2] * C_1_tp + u_tp_fp[,2]
lambda_2_fp <- 1 / (1 + exp(-tmp))
C_2_tp <- rbinom(N, 1, prob = lambda_2_tp)
C_2_fp <- rbinom(N, 1, prob = lambda_2_fp)
for(i in 1:N){
  pi_tp_tmp <- alpha_tp[1,] + C_2_tp[i] * alpha_tp[2,]
  y_tp[N+i,] <- rbinom(J, 1, prob = 1 / (1 + exp(-pi_tp_tmp)))
  pi_fp_tmp <- alpha_fp[1,] + C_2_fp[i] * alpha_fp[2,]
  y_fp[N+i,] <- rbinom(J, 1, prob = 1 / (1+ exp(-pi_fp_tmp)))
}
## generate responses for wave 3
tmp <- xcovs %*% b_tp[,3] + beta_tp[1] * C_2_tp + beta_tp[2] * C_2_fp + u_tp_fp[,1]
lambda_3_tp <- 1 / (1 + exp(-tmp))
tmp <- xcovs %*% b_fp[,3] + beta_fp[1] * C_2_fp + beta_fp[2] * C_2_tp + u_tp_fp[,2]
lambda_3_fp <- 1 / (1 + exp(-tmp))
C_3_tp <- rbinom(N, 1, prob = lambda_3_tp)
C_3_fp <- rbinom(N, 1, prob = lambda_3_fp)
for(i in 1:N){
  pi_tp_tmp <- alpha_tp[1,] + C_3_tp[i] * alpha_tp[2,]
  y_tp[2*N+i,] <- rbinom(J, 1, prob = 1 / (1+exp(-pi_tp_tmp)))
  pi_fp_tmp <- alpha_fp[1,] + C_3_fp[i] * alpha_fp[2,]
  y_fp[2*N+i,] <- rbinom(J, 1, prob = 1 / (1+exp(-pi_fp_tmp)))
}
C_it_tp <- c(C_1_tp, C_2_tp, C_3_tp)
C_it_fp <- c(C_1_fp, C_2_fp, C_3_fp)
colMeans(matrix(C_it_tp,N))
colMeans(matrix(C_it_fp,N))

alpha_tp0 <- alpha_fp0 <- rbind(rep(0, J), rep(0.5, J))
struc_params0 <- list(b_tp = b_tp,
                      b_fp = b_fp,
                      d_tp = d_tp,
                      d_fp = d_fp,
                      beta_tp = beta_tp,
                      beta_fp = beta_fp,
                      rho_u = rho_u)
res <- latent_trans2p(y_tp, y_fp, C_it_tp, C_it_fp, u_tp_fp[,1], u_tp_fp[,2],
                      alpha_tp0, alpha_fp0, struc_params0, mcmc_len = 1000)
alpha_tp
matrix(colMeans(res$alpha_tp_draws),2)
alpha_fp
matrix(colMeans(res$alpha_fp_draws),2)
cbind(unlist(struc_params0),colMeans(res$struc_params_draws))
plot(res$struc_params_draws[,13])
plot(res[[1]][,1])


struc_paramdraws <- matrix(0, 100, length(unlist(struc_params0)))
for(i in 1:100){
  cat("i= ", i, "\n")
  struc_params0 <- sample_struc2p(C_it_tp, C_it_fp, struc_params0)
  struc_paramdraws[i,] <- unlist(struc_params0)
}
cbind(unlist(struc_params0),colMeans(struc_paramdraws))
# estimated measurement parameters
matrix(colMeans(res[[1]]),2)
# true measurement parameters
alpha_tp

# estimated & true structural parameters
struc_est <- rbind(colMeans(res[[2]]),c(b_tp, d_tp, beta_tp, sig2_u_tp))
rownames(struc_est) <- c('estimated', 'true')
colnames(struc_est) <- c('b1','b2','b3', 'd_tp', 'beta', 'sigma2')
struc_est

plot(1:1000, res[[2]][,6], 'l')

plot(1:1000, (1:1000)^2,'l')
