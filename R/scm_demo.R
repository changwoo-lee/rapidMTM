#### Demo of Bayesian SCM       ####
#### With multiple-try proposal ####

rm(list = ls())
source("scm_mtm.R")
source("scm_weight_functions.R")

### Generate data -----

# read in piecewise constants
raw_data = R.matlab::readMat("cluster_scm.mat")
coords = cbind(raw_data$lon, raw_data$lat)
colnames(coords) = c("lon", "lat")
mu_true = raw_data$beta[, 2]
cluster_true = as.integer(as.factor(mu_true))
n = length(mu_true)

# get minimum spanning tree from Delaunay triangulation graph
graph0 = dentrigraph(coords, threshold = 0.1)
mstgraph0 = mst(graph0)
if('weight' %in% names(edge_attr(mstgraph0)))
  mstgraph0 = delete_edge_attr(mstgraph0, 'weight')

# look at the true cluster and graph
plotGraph(coords, mstgraph0, title = "Spanning tree graph and true clusters") +
  geom_point(aes(x = lon, y = lat, color = factor(cluster_true)), 
             size = 1, data = as.data.frame(coords))

# check how many edges should be removed from mstgraph0 to get the true clusters
estatus_mst = getEdgeStatus(cluster_true, as_edgelist(mstgraph0, names = F))
eid_btw_mst_true = which(estatus_mst == 'b')
length(eid_btw_mst_true)

# get true compatible partition
mst_subgraph = delete.edges(mstgraph0, eid_btw_mst_true)
cluster_true_comp = components(mst_subgraph)$membership

# generate y
set.seed(1234)
snr = 10
Y = mu_true + rnorm(n, 0, sd(mu_true) / snr)

# hyperparameters
hyperpar = c()
hyperpar["a0"] = 1
hyperpar["b0"] = 1
hyperpar["c0"] = 0.5
hyperpar["lambda"] = 0.01

# compute posterior of the true clustering state (compatible to spanning tree)
csize_true = as.integer(table(cluster_true_comp))
csum_true = as.numeric(tapply(Y, cluster_true_comp, sum))
log_post_true = evalLogPost(n, sum(Y ^ 2), csize_true, csum_true, hyperpar)


### Fit BSCM with multiple try -----

## set weight function (see scm_weight_functions.R)
# lweight_fcn = lweight_fcn_sqrt
# lweight_fcn = lweight_fcn_ord
lweight_fcn = lweight_fcn_min
# lweight_fcn = lweight_fcn_max

## set number of trails
MT = 100

## initialize with SCC
mu_scc = fitSCC(Y, mstgraph0)
init_cluster = as.integer(as.factor(mu_scc))

MCMC = 1000
BURNIN = 0
THIN = 1

## fit BSCC-MT (see scm_mtm.R for detailed documentation)
mcmc_res = fitBSCM(Y, mstgraph0, init_cluster, hyperpar, lweight_fcn, MCMC, BURNIN, THIN, MT, 
                   seed = 1234, cluster_true = cluster_true_comp)

## traceplot of log posterior
plot(
  mcmc_res$log_post_out, 
  type = "l", 
  main = "Trace plot of log posterior", 
  xlab = "iteration", 
  ylab = "log posterior (up to a constant)"
)
# add true state log posterior
abline(h = log_post_true, col = "red", lty = 2)
