library(armspp)
library(MCMCpack)
library(hmmgibbs)

### simulate single process data
N <- 16439
wave_num <- 3
J <- 8

## coefficients for pi
alpha_tp <- rbind(rnorm(J, 0,1), runif(J))

## random effect u_tp, u_fp
sig2_u_tp <- 1
u_tp <- rnorm(N, 0, sd = sqrt(sig2_u_tp))

## time variant intercepts for lambda
b_tp <- c(-1, 0, 1)

## coefficients for cross-lagged effect
beta_tp <- 0.5

## scale parameter for random effects in wave 1
d_tp <- 0.5


y_tp <- matrix(0, N*wave_num, J)
tind <- rep(1:3,N)
ind <- rep(1:N, each=3)

## generate responses for wave 1
## generate lambda_1^{tp}
tmp <- b_tp[1] + d_tp * u_tp
lambda_1_tp <- 1 / (1+exp(-tmp))

## generate latent class variables C_1_tp, C_1_fp
C_1_tp <- rbinom(N, 1, prob = lambda_1_tp)

for(i in 1:N){
  pi_tp_tmp <- alpha_tp[1,] + C_1_tp[i] * alpha_tp[2,]
  y_tp[i,] <- rbinom(J, 1, prob = 1 / (1 + exp(-pi_tp_tmp)))
}
## generate responses for wave 2
tmp <- b_tp[2] + beta_tp[1] * C_1_tp + u_tp
lambda_2_tp <- 1 / (1 + exp(-tmp))
C_2_tp <- rbinom(N, 1, prob = lambda_2_tp)
for(i in 1:N){
  pi_tp_tmp <- alpha_tp[1,] + C_2_tp[i] * alpha_tp[2,]
  y_tp[N+i,] <- rbinom(J, 1, prob = 1 / (1 + exp(-pi_tp_tmp)))
}
## generate responses for wave 3
tmp <- b_tp[3] + beta_tp[1] * C_2_tp + u_tp
lambda_3_tp <- 1 / (1 + exp(-tmp))
C_3_tp <- rbinom(N, 1, prob = lambda_3_tp)
for(i in 1:N){
  pi_tp_tmp <- alpha_tp[1,] + C_3_tp[i] * alpha_tp[2,]
  y_tp[2*N+i,] <- rbinom(J, 1, prob = 1 / (1+exp(-pi_tp_tmp)))
}
C_it_tp <- c(C_1_tp, C_2_tp, C_3_tp)
C_mat <- matrix(C_it_tp, N)
colMeans(C_mat)
# res <- sample_struc(C_it_tp, u_tp, c(rep(1,3), 1, 1, 1), 1000)
# plot(1:1000, res[,5], 'l')
# sample b1-bt, d1, beta, sigma2

alpha_tp0 <- rbind(rep(0, J), rep(0.5, J))
# b1-bt, d1, beta, sigma2
struc_params0 <- c(rep(0,3), 1, 1, 1)
free_params <- c(1,1,1,1,1,1)

res <- latent_trans(y_tp, C_it_tp, u_tp,
                    alpha_tp0, struc_params0, free_params, 1000)

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
