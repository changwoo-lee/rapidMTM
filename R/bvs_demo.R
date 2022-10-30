# bvs_demo.R
rm(list=ls())
## Load libraries, if not installed, install with: install.packages("package name")
library(mgcv) # for Cholesky update
library(matrixStats) # for logsumexp

## data generation
n = 1000; p = 5000; # dimension of design matrix
snr = 4 # signal-to-noise ratio

s = 10; b0 = c(2,-3,2,2,-3,3,-2,3,-2,3); # true model size = 10
beta = c(b0, rep(0, p - s)) * sqrt( log(p)/n ) * snr # true coefficients

set.seed(1)
X = scale(matrix(rnorm(n*p), nrow=n), center = T)
y = X %*% beta + rnorm(n)

# preprocess, take some time
preprocessed = list(
  XtX = crossprod(X),
  Xty = crossprod(X, y),
  yty = crossprod(y)
)

# initialize gamma (binary vector) s.t. 20 steps needed to reach the true model 
gammainit = rep(0,p) 
gammainit[sample.int(p-10, size = 10) + 10] = 1

# under hyperparameter settings kappa = 2, g = p^3
kappa = 2; g = p^3; 
strue = sum(beta!=0)
SSRtrue = preprocessed$yty - g/(g+1)*t(preprocessed$Xty[1:10])%*%solve(preprocessed$XtX[1:10,1:10],preprocessed$Xty[1:10])
logposttrue = -kappa*strue*log(p) - strue/2*log(1+g) - n/2*log(SSRtrue) # see eq. 13 of the paper
# run MTM 

source("R/bvs_mtm.R") # MTM algorithm for BVS
# ntry = 20  
fit_ord_20 = bvs_mtm(y, X, ntry =  20, balancingft = NULL, burn = 0, nmc = 2000,
                     preprocessed = preprocessed, gammainit = gammainit)

fit_sqrt_20 = bvs_mtm(y, X, ntry =  20, balancingft = "sqrt", burn = 0, nmc = 2000,
                      preprocessed = preprocessed, gammainit = gammainit)

fit_min_20 = bvs_mtm(y, X, ntry =  20, balancingft = "min", burn = 0, nmc = 2000,
                     preprocessed = preprocessed, gammainit = gammainit)

fit_max_20 = bvs_mtm(y, X, ntry =  20, balancingft = "max", burn = 0, nmc = 2000,
                     preprocessed = preprocessed, gammainit = gammainit)

# ntry = 200
fit_ord_200 = bvs_mtm(y, X, ntry =  200, balancingft = NULL, burn = 0, nmc = 2000,
                      preprocessed = preprocessed, gammainit = gammainit)

fit_sqrt_200 = bvs_mtm(y, X, ntry =  200, balancingft = "sqrt", burn = 0, nmc = 2000,
                       preprocessed = preprocessed, gammainit = gammainit)

fit_min_200 = bvs_mtm(y, X, ntry =  200, balancingft = "min", burn = 0, nmc = 2000,
                      preprocessed = preprocessed, gammainit = gammainit)

fit_max_200 = bvs_mtm(y, X, ntry =  200, balancingft = "max", burn = 0, nmc = 2000,
                      preprocessed = preprocessed, gammainit = gammainit)

# plot 
plot(fit_ord_20$logpostout, type="l", main = "Trace plot of log posterior", xlab = "iteration", ylab="log posterior (up to a constant)")
lines(fit_sqrt_20$logpostout, col = 2)
lines(fit_min_20$logpostout, col = 3)
lines(fit_max_20$logpostout, col = 7)
lines(fit_ord_200$logpostout, col = 1, lty = 2)
lines(fit_sqrt_200$logpostout, col = 2, lty = 2)
lines(fit_min_200$logpostout, col = 3, lty = 2)
lines(fit_max_200$logpostout, col = 7, lty = 2)
abline(h = logposttrue, col = 4, lty = 2, lwd =3)
legend("bottomright",lty = c(1,1,1,1,2,2,2,2,2),
       col = c(1,2,3,7,1,2,3,7,4),lwd = c(1,1,1,1,1,1,1,1,3),
       legend = c("ord, N=20","sqrt, N=20","min, N=20","max, N=20",
                  "ord, N=200","sqrt, N=200","min, N=200","max, N=200",
                  "true model"))

### ntry = 2000
fit_ord_2000 = bvs_mtm(y, X, ntry =  2000, balancingft = NULL, burn = 0, nmc = 1000,
                       preprocessed = preprocessed, gammainit = gammainit)

fit_sqrt_2000 = bvs_mtm(y, X, ntry =  2000, balancingft = "sqrt", burn = 0, nmc = 1000,
                        preprocessed = preprocessed, gammainit = gammainit)

fit_min_2000 = bvs_mtm(y, X, ntry =  2000, balancingft = "min", burn = 0, nmc = 1000,
                       preprocessed = preprocessed, gammainit = gammainit)

fit_max_2000 = bvs_mtm(y, X, ntry =  2000, balancingft = "max", burn = 0, nmc = 1000,
                       preprocessed = preprocessed, gammainit = gammainit)

plot(fit_sqrt_2000$logpostout, type = "l", main = "Trace plot of log posterior", xlab = "iteration", ylab="log posterior (up to a constant)")
lines(fit_ord_2000$logpostout, col = 1)
lines(fit_min_2000$logpostout, col = 3)
lines(fit_max_2000$logpostout, col = 7)
abline(h = logposttrue, col = 4, lty = 2, lwd =3)
legend("right", lty = c(1,1,1,1,2), col = c(1,2,3,7,4),lwd = c(1,1,1,1,3),
       legend = c("ord, N=2000","sqrt, N=2000","min, N=2000","max, N=2000",
                  "true model"))


