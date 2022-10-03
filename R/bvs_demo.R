# bvs_demo.R
# Part of the supplementary material:
# Anonymous authors, Rapidly Mixing Multiple-try Metropolis Algorithms for Model Selection Problems.
# submitted to NeurIPS 2022, paper7322

source("bvs_singletry.R")
source("bvs_mtm_ord.R")
source("bvs_mtm_sqrt.R")
source("bvs_mtm_min.R")
source("bvs_mtm_max.R")

snr = 4
n = 1000
p = 5000

s = 10
b0 = c(2,-3,2,2,-3,3,-2,3,-2,3)
beta = c(b0, rep(0, p - s)) * sqrt( log(p)/n ) * snr # noise sd = 1

set.seed(1)

# if indep
X = matrix(rnorm(n*p), nrow=n)
# if corr
# 
# C = diag(p)
# for (i in 1:p){
#   for (j in 1:p){
#     C[i,j] = exp(-abs(i - j))
#   }
# }
# R = chol(C)
# 
# X = matrix(rnorm(n*p), nrow=n) %*% R

X = scale(X, center=TRUE)
y = X %*% beta + rnorm(n)

# preprocess, take some time
preprocessed = list(
  XtX = crossprod(X),
  Xty = crossprod(X, y),
  yty = crossprod(y)
)

# initial gamma (variables)
gammainit = rep(0,p)
gammainit[sample.int(p-10, size = 10) + 10] = 1

# fit single-try MH
fit_single = bvs_singletry(y, X, g = p^3, kappa = 2, s0= 100, 
                           burn = 0, nmc = 1e5, thin = 1,
                           gammainit = gammainit, 
                           verbose = F, debug = F, preprocessed = preprocessed, 
                           truegamma = as.logical(beta)) # measuring 'hit'

fit_ordinary <- bvs_mtm_ord(y = y, X = X, g = p^3, kappa = 2, s0 = 100,
                        burn = 0, nmc = 1e3, thin = 1,
                        K= 100, ########### number of tries
                        gammainit = gammainit, preprocessed = preprocessed, 
                        truegamma = as.logical(beta)) # measuring 'hit'

fit_sqrt <- bvs_mtm_sqrt(y = y, X = X, g = p^3, kappa = 2, s0 = 100,
                                 burn = 0, nmc = 1e3, thin = 1,
                                 K= 100, ########### number of tries
                                 gammainit = gammainit, preprocessed = preprocessed, 
                                 truegamma = as.logical(beta)) # measuring 'hit'


# true log posterior prob.
g = p^3
kappa = 2
SSRtrue = preprocessed$yty - g/(g+1)*preprocessed$Xty[1:s]%*%solve(preprocessed$XtX[1:s,1:s],preprocessed$Xty[1:s])
truelogpost = lgamma(n/2)- (n/2)*log(pi) - s/2*log(1+g) - n/2*log(SSRtrue) - kappa*s*log(p)
truelogpost = as.numeric(truelogpost)


plot(fit_single$logpostout, type="l", main = "single-try MH")
abline(h=truelogpost, col = 2)

plot(fit_ordinary$logpostout, type="l", main = "ordinary, ntry=100")
abline(h=truelogpost, col = 2)

plot(fit_sqrt$logpostout, type="l", main = "sqrt, ntry=100")
abline(h=truelogpost, col = 2)


## if ntry is too large, e.g. ntry = 2000
fit_ordinary_2000 <- bvs_mtm_ord(y = y, X = X, g = p^3, kappa = 2, s0 = 100,
                            burn = 0, nmc = 1e3, thin = 1,
                            K= 2000, ########### number of tries
                            gammainit = gammainit, preprocessed = preprocessed, 
                            truegamma = as.logical(beta)) # measuring 'hit'

fit_sqrt_2000 <- bvs_mtm_sqrt(y = y, X = X, g = p^3, kappa = 2, s0 = 100,
                                 burn = 0, nmc = 1e3, thin = 1,
                                 K= 2000, ########### number of tries
                                 gammainit = gammainit, preprocessed = preprocessed, 
                                 truegamma = as.logical(beta)) # measuring 'hit'

plot(fit_sqrt_2000$logpostout, type="l", main = "sqrt(solid) and ordinary(dash), ntry=2000")
lines(fit_ordinary_2000$logpostout, lty=2) # dash
abline(h=truelogpost, col = 2)

