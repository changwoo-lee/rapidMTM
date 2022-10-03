# sbm_demo.R 
# Part of the supplementary material:
# Anonymous authors, Rapidly Mixing Multiple-try Metropolis Algorithms for Model Selection Problems.
# submitted to NeurIPS 2022, paper7322

source("sbm_singletry.R")
source("sbm_mtm_sqrt.R")
source("sbm_mtm_ord.R")
source("sbm_mtm_min.R")
source("sbm_mtm_max.R")

library(igraph)

# weak: CH = 2.0
n= 1000
k= 2
a = 0.07  # within-prob
b = 0.01 # cross-prob
n/(k*log(n))*(sqrt(a)-sqrt(b))^2 # around 2



set.seed(1)

Q = matrix(b, k, k)
diag(Q) = a
g = sample_sbm(n, pref.matrix = Q, block.sizes = rep(n/k, k))
components(g)$no # should be one component
true_z = rep(1:k, each = n/k)
A = as_adj(g, sparse = F)


#initial partition z such that hamming distance between the true is 400
z_init = true_z
wrongidx = sample.int(n, size = 2*n/5)
z_init[wrongidx] = (true_z[wrongidx] + sample(k-1, size = 2*n/5, replace = T))%% k
z_init[which(z_init==0)] = k

# to check d_H = 400 is correct, please uncomment below and run
# library(partitionComparison)
# n*(partitionComparison::classificationErrorDistance(new("Partition", z_init),new("Partition", true_z)))


fit_single = sbm_singletry(A, k, z_init = z_init, 
                           N_iter = 200000, true_z = true_z)

fit_ordinary = sbm_mtm_ord(A, k, ntry = 100, 
                           z_init = z_init, 
                           N_iter = 2000, 
                           true_z = true_z, sparse = F)

fit_sqrt = sbm_mtm_sqrt(A, k, ntry = 100, 
                        z_init = z_init, 
                        N_iter = 2000, 
                        true_z = true_z, sparse = F)


truelogpost = calculatelog_lik(A, true_z)



plot(fit_single$loglik_post, type="l", main = "single-try MH")
abline(h=truelogpost, col = 2)

plot(fit_sqrt$loglik_post, type="l", main = "sqrt(solid) and ordinary(dash), ntry=100")
lines(fit_ordinary$loglik_post, lty=2) # dash
abline(h=truelogpost, col = 2)




