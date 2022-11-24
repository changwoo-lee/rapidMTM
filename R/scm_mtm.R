### Spatial clustering model using finxed spanning tree graph  ###
### With multiple-try proposal                                 ###

require(igraph)
library(deldir)
require(genlasso) 
require(ggplot2)
require(glmnet)

### Functions for fitting SCM with MTM -----

#' fitBSCC is the main function to run MTM with a custom weight function for Bayesian spatial clustering problems
#'
#' @param Y Vector of responses of length n
#' @param stgraph Spanning tree (igraph object) connecting n spatial locations
#' @param init_cluster Initial cluster assignments of MCMC. Should be an n-dimensional vector where the i-th element is the id of the cluster that the i-th location belongs to.
#'                     Each cluster id should be an integer ranging from 1 to the number of clusters.
#' @param hyper Hyperparameters of MCMC. Should be a named vector with the following elements:
#'              "a0", "b0", "c0", "lambda": Prior parameters
#'              "k_max": (Optional) maximum umber of clusters
#' @param lweight_fcn Weight function returning logarithm of weights. Should be a function that accepts the following two parameters:
#'                    The first parameter is the log posterior of a proposed state.
#'                    The second parameter is the log posterior of the current state.
#'                    See \code{scm_weight_functions.R} for example weight functions.
#' @param MCMC Number of MCMC iterations. The number of poseterior samples drawn will be (MCMC-BURNIN)/THIN
#' @param BURNIN Number of burn-in iterations
#' @param THIN Size of thinning interval
#' @param MT Number of trials in MTM. Default is 1 (Metropolis-Hastings)
#' @param seed Random seed. Default is 1234
#' @param cluster_true (Optional) n-dimensional vector of true cluster memberships used to compute hitting time and iterations. 
#'
#' @return A list of the following elements, where n_post is the number of posterior samples:
#'         "cluster_out": n_post * n matrix of cluster memberships after burn-in and thinning
#'         "log_post_out": n_post-dimensional vector of log posterior after burn-in and thinning
#'         "log_post_init": log posterior of initial state
#'         "cnt_acc": Number of accepted proposals
#'         "time_mcmc": Wall time of the MCMC run
#'         If \code{cluster_true} is provided:
#'         "rand_out": n_post-dimensional vector of Rand index of each posterior sample
#'         "iter_hit": First iteration that hits the true clustering state (ignoring burn-in and thinning). Contain Inf if the truth is never hit.
#'         "time_hit": Wall time to hit the true clustering state. Contain Inf if the truth is never hit.
#'         "iter_hit_099": First iteration that hits the 0.99 Rand index neighbor of the true clustering state (ignoring burn-in and thinning). Contain Inf if the truth is never hit.
#'         "time_hit_099": Wall time to hit the 0.99 Rand index neighbor of the true clustering state. Contain Inf if the truth is never hit.
fitBSCM <- function(Y, stgraph, init_cluster, hyper, lweight_fcn, MCMC, BURNIN, THIN, MT = 1, seed = 1234,
                    cluster_true = NULL) {
  set.seed(seed)
  
  n = vcount(stgraph)
  YtY = sum(Y ^ 2)
  # hyper-parameter
  k_max = ifelse(is.na(hyper['k_max']), n, hyper['k_max'])
  
  if('name' %in% names(vertex_attr(stgraph))) {
    stgraph = delete_vertex_attr(stgraph, 'name')
  }
  V(stgraph)$vid = 1:vcount(stgraph)
  E(stgraph)$eid = 1:ecount(stgraph)
  inc_mat = get.edgelist(stgraph, names = F)
  
  # initial values
  cluster = init_cluster  # vector of length n
  k = max(cluster)  # number of clusters
  
  csize = as.integer(table(cluster))  # cluster size
  csum = as.numeric(tapply(Y, cluster, sum))  # sum of Y within each cluster
  
  # get induced subgraph of stgraph for each cluster
  clust_vid = split(1:n, cluster)
  subgraphs = lapply(
    clust_vid, 
    function(vids, stgraph) induced_subgraph(stgraph, vids), 
    stgraph
  )
  
  # index for between-cluster edges
  c1 = cluster[inc_mat[, 1]]; c2 = cluster[inc_mat[, 2]]
  idx_btw = which(c1 != c2)
  eid_btw = (E(stgraph)$eid)[idx_btw]
  
  # log marginal posterior
  log_post = evalLogPost(n, YtY, csize, csum, hyper)
  log_post_init = log_post
  
  ## RUn MCMC ##
  
  ## MCMC results
  cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, n))
  log_post_out = numeric((MCMC-BURNIN)/THIN)
  rand_out = numeric((MCMC-BURNIN/THIN))
  
  cnt_acc = 0
  
  # timing
  time_mcmc = Inf
  time_hit = Inf
  iter_hit = Inf
  hit = F
  time_hit_099 = Inf    # time to reach rand index 0.99
  iter_hit_099 = Inf    # iteration to reach rand index 0.99
  hit_099 = F           # hit rand index 0.99?
  time_start = Sys.time()
  
  ## MCMC iteration
  for(iter in 1:MCMC) {
    proposal_mt = lapply(
      1:MT,
      FUN = function(i) proposePartition(k, Y, YtY, hyper, log_post, lweight_fcn,
                                         subgraphs, csize, csum, cluster, eid_btw, inc_mat)
    )
    
    # sample one proposal based on (log) weights
    lweights = sapply(proposal_mt, FUN = function(res) res$lweight)
    lweights_normalized = lweights - logSumExp(lweights)
    proposal = proposal_mt[[sample.int(MT, 1, prob = exp(lweights_normalized))]]
    log_post_new = proposal$log_post
    csum_new = proposal$csum
    move = proposal$move
    rm(proposal_mt)
    
    if(move == 'split') { ## Birth move
      
      # update information for the selected proposal
      update_res = updateSplit(proposal, subgraphs, eid_btw)
      subgraphs_new = update_res$subgraphs
      csize_new = proposal$csize
      eid_btw_new = update_res$eid_btw
      cluster_new = update_res$cluster
      k_new = k + 1
      
    } else { ## Death move
      
      # update information for the selected proposal
      update_res = updateMerge(proposal, subgraphs, eid_btw, stgraph)
      subgraphs_new = update_res$subgraphs
      csize_new = proposal$csize
      eid_btw_new = update_res$eid_btw
      cluster_new = update_res$cluster
      k_new = k - 1
      
    }
    
    proposal_mt_new = lapply(
      1:(MT-1),
      FUN = function(i) proposePartition(k_new, Y, YtY, hyper, log_post_new, lweight_fcn,
                                         subgraphs_new, csize_new, csum_new, cluster_new, eid_btw_new, inc_mat)
    )
    
    lweights_new = sapply(proposal_mt_new, FUN = function(res) res$lweight)
    lweights_new[MT] = lweight_fcn(log_post, log_post_new)
    rm(proposal_mt_new)
    
    # acceptance probability
    acc_prob = min(0, logSumExp(lweights) - logSumExp(lweights_new))
    acc_prob = exp(acc_prob)
    if(runif(1) < acc_prob){
      # accept
      subgraphs = subgraphs_new
      csize = csize_new
      csum = csum_new
      eid_btw = eid_btw_new
      cluster = cluster_new
      k = k_new
      log_post = log_post_new
      
      cnt_acc = cnt_acc + 1
    }
    
    
    ## save result
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      cluster_out[(iter-BURNIN)/THIN, ] = cluster
      log_post_out[(iter-BURNIN)/THIN] = log_post
      
      # hitting the truth?
      if (!is.null(cluster_true)) {
        rand = fossil::rand.index(cluster_true, cluster)
        rand_out[[(iter-BURNIN)/THIN]] = rand
        if (!hit & rand == 1) {
          # hit!
          iter_hit = iter
          time_hit = difftime(Sys.time(), time_start, units = "secs")
          hit = T
        }
        
        if (!hit_099 & rand > 0.99) {
          # hit!
          iter_hit_099 = iter
          time_hit_099 = difftime(Sys.time(), time_start, units = "secs")
          hit_099 = T
        }
      }
    }
    
    if(iter %% 100 == 0) cat('Iteration', iter, 'done\n')
  }
  time_mcmc = difftime(Sys.time(), time_start, units = "secs")
  
  mode(cluster_out) = 'integer'  # to save memory
  return(list('cluster_out' = cluster_out, 'log_post_out' = log_post_out, 'log_post_init' = log_post_init,
              'cnt_acc' = cnt_acc, 'rand_out' = rand_out,
              'iter_hit' = iter_hit, 'time_hit' = time_hit, 'time_mcmc' = time_mcmc,
              'iter_hit_099' = iter_hit_099, 'time_hit_099' = time_hit_099))
}


# function to split an existing cluster given a spanning tree
# subgraphs: cluster_id -> subgraph
# cluster: vid -> cluster_id
splitCluster <- function(subgraphs, csize, cluster, clust_split = NULL, edge_cutted = NULL) { 
  k = length(subgraphs)
  
  if (is.null(clust_split) | is.null(edge_cutted)) {
    clust_split = sample.int(k, 1, prob = csize - 1)
    st_subgraph = subgraphs[[clust_split]]
    
    eidx = sample.int(csize[clust_split]-1, 1)
    edge_cutted = E(st_subgraph)[eidx]
  } else {
    st_subgraph = subgraphs[[clust_split]]
  }
  
  eid_cutted = edge_cutted$eid
  st_subgraph = delete.edges(st_subgraph, edge_cutted)
  connect_comp = components(st_subgraph)
  idx_new = (connect_comp$membership == 2)
  vid_new = V(st_subgraph)$vid[idx_new]
  vid_old = V(st_subgraph)$vid[!idx_new]
  
  cluster[vid_new] = k + 1
  csize[clust_split] = length(vid_old)
  csize[k + 1] = length(vid_new)
  
  return(list(vid_old = vid_old, vid_new = vid_new, eid_cutted = eid_cutted, csize = csize,
              clust_old = clust_split, idx_new = idx_new, cluster_new = cluster, k = k + 1))
}

# function to update if a split move is accepted
updateSplit <- function(split_res, subgraphs, eid_btw) {
  k = length(subgraphs)
  
  clust_split = split_res$clust_old
  vid_old = split_res$vid_old; vid_new = split_res$vid_new
  
  subgraph_split = subgraphs[[clust_split]]
  idx_new = split_res$idx_new
  subgraphs[[clust_split]] = induced_subgraph(subgraph_split, !idx_new)  # subgraph of old cluster
  subgraphs[[k+1]] = induced_subgraph(subgraph_split, idx_new) # subgraph of new cluster
  
  cluster = split_res$cluster
  eid_btw = c(eid_btw, split_res$eid_cutted)
  
  return(list(subgraphs = subgraphs, cluster = cluster, eid_btw = eid_btw))
}

# function to merge two existing clusters
mergeCluster <- function(eid_btw, subgraphs, csize, cluster, edge_list, eidx = NULL) {
  # edge for merging
  edge_merge = ifelse(is.null(eidx), sample.int(length(eid_btw), 1), eidx)
  # update cluster information
  # clusters of endpoints of edge_merge
  eid_merge = eid_btw[edge_merge]
  clusters_merge = cluster[edge_list[eid_merge, ]]
  clusters_merge = sort(clusters_merge)
  c1 = clusters_merge[1]; c2 = clusters_merge[2] # note c1 < c2
  # merge c2 to c1
  
  # vid of vertices in c2
  vid_old = V(subgraphs[[c2]])$vid
  # vid in merged cluster
  vid_new = c(V(subgraphs[[c1]])$vid, vid_old)
  
  # update cluster memberships
  cluster[vid_old] = c1
  idx = which(cluster > c2)
  cluster[idx] = cluster[idx] - 1
  
  # update csize
  csize[c1] = length(vid_new)
  csize = csize[-c2]
  
  # now drop c2
  return(list(vid_old = vid_old, vid_new = vid_new, clust_old = c2, clust_new = c1, csize = csize,
              edge_merge = edge_merge, cluster = cluster, k = length(subgraphs) - 1))
}

# function to update if a merge move is accepted
updateMerge <- function(res_merge, subgraphs, eid_btw, stgraph) {
  clust_old = res_merge$clust_old; clust_new = res_merge$clust_new
  vid_old = V(subgraphs[[clust_old]])$vid
  vid_new = c(V(subgraphs[[clust_new]])$vid, vid_old)
  subgraphs[[clust_new]] = induced_subgraph(stgraph, vid_new)
  subgraphs[[clust_old]] = NULL
  
  cluster = res_merge$cluster
  
  eid_btw = eid_btw[-res_merge$edge_merge]
  
  return(list(subgraphs = subgraphs, cluster = cluster, eid_btw = eid_btw))
}

# function to propose a move to modify partition
proposePartition <- function(k, Y, YtY, hyper, log_post, lweight_fcn, 
                             subgraphs, csize, csum, cluster, eid_btw, edge_list) {
  n = length(Y)
  split_prob = (n - k) / (n - 1)
  
  if (runif(1) < split_prob) {
    move = 'split'
    proposal = splitCluster(subgraphs, csize, cluster)
    
    # update csum
    csum_old = csum[proposal$clust_old]
    csum[proposal$clust_old] = sum(Y[proposal$vid_old])
    csum[k + 1] = csum_old - csum[proposal$clust_old]
    
  } else {
    move = 'merge'
    proposal = mergeCluster(eid_btw, subgraphs, csize, cluster, edge_list)
    
    # update csum
    csum[proposal$clust_new] = csum[proposal$clust_new] + csum[proposal$clust_old]
    csum = csum[-proposal$clust_old]
  }
  
  # compute log posterior
  log_post_new = evalLogPost(n, YtY, proposal$csize, csum, hyper)
  proposal[['log_post']] = log_post_new
  
  # compute weights
  lweight = lweight_fcn(log_post_new, log_post)
  proposal[['lweight']] = lweight
  
  proposal[['move']] = move
  proposal[['csum']] = csum
  
  return(proposal)
}

# function to get log marginal t likelihood (up to a constant)
evalLogMarLike <- function(n, ete, csize, csum, lambda, a0, b0) {
  k = length(csize)
  Xte = csum
  quad = ete - sum(Xte^2 / (lambda + csize))
  log_like = -(n + a0) / 2 * log((b0 + quad) / 2)
  log_det = -(1/2) * (sum(log(lambda + csize)) - k * log(lambda))
  log_like = log_like + log_det
  return(log_like)
}


# function to get log marginal posterior density (up to a constant)
evalLogPost <- function(n, YtY, csize, csum, hyper) {
  k = length(csize)
  a0 = hyper['a0']; b0 = hyper['b0']; c0 = hyper['c0']
  lambda = hyper['lambda']
  log_prior = -lchoose(n-1, k-1) + k * log(1 - c0)
  log_like = evalLogMarLike(n, YtY, csize, csum, lambda, a0, b0)
  return(log_like + log_prior)
}

### Funtions for fitting SCC (Li and Sang, 2019, JASA) -----

calInvDMST <- function(n, mstgraph) {
  edgelist = get.edgelist(mstgraph)
  D = getDgSparse(mstgraph)
  D = rbind2(D, rep(1/n, n))
  inv_D = solve(D)
  return(inv_D)
}

convertX <- function(inv_D, Xi, n, p) {
  X.mat.trans = matrix(0, n, n*p)
  for(j in 1:p) { 
    X.mat.trans[, ((j-1)*n+1):(j*n)] = as.matrix(Xi[, j] * inv_D);
  }
  return(X.mat.trans)
}

convertBeta <- function(inv_D, theta, n, p){
  beta_hat = matrix(0, nrow(theta), ncol(theta))
  for(i in 1:p) {
    beta_hat[((i-1)*n+1):(i*n), ] = as.matrix(inv_D %*% theta[((i-1)*n+1):(i*n), ])
  }
  return(matrix(beta_hat, n, p))
}

estBetaBIC <- function(Y, X, int.index, p, weight, lambda, inv_D) {
  path.fit.0 = glmnet(X[, -int.index], Y, family = "gaussian", standardize = FALSE, lambda = lambda, 
                      intercept = TRUE, penalty.factor = weight[-int.index])
  path.fit = path.fit.0
  path.fit$beta = rbind(path.fit.0$beta, path.fit.0$a0)
  path.fit$beta[-int.index, ] = path.fit.0$beta
  path.fit$beta[int.index, ] = path.fit.0$a0
  
  IC.id = selLambdaBIC(Y, X, path.fit)

  est.Beta = list(
    "BIC" = convertBeta(inv_D, path.fit$beta[, IC.id$BIC, drop = FALSE], n, p)
  )
  return(est.Beta)
}

selLambdaBIC <- function(Y, X, fit.lasso) {
  n = length(Y)
  p = ncol(X) / n;
  k = fit.lasso$df;
  fit.beta = fit.lasso$beta;
  MSE = apply((Y - X %*% fit.beta) ^ 2, 2, mean) 
  BIC = n * log(MSE) + k * log(n); ## plug in estimates of sigma^2
  return(list("BIC" = which.min(BIC)))
}

# function to fit SCC (Li and Sang, 2019, JASA)
fitSCC <- function(Y, mstgraph0) {
  n = length(Y)
  inv_D = calInvDMST(n, mstgraph0)
  X.mat.trans = convertX(inv_D, matrix(rep(1, n), ncol = 1), n, 1)
  weight = rep(c(rep(1, n-1), 0), 1) 
  lambda = 10 ^ seq(-6, -2, length = 300) 
  est.Beta.IC = estBetaBIC(Y, X.mat.trans, n, 1, weight, lambda, inv_D)
  mu_scc = est.Beta.IC$BIC[, 1]
  return(mu_scc)
}

### Utils -----

# function to get whether an edge is within a cluster or between two clusters
getEdgeStatus <- function(membership, inc_mat) {
  membership_head = membership[inc_mat[, 1]]
  membership_tail = membership[inc_mat[, 2]]
  edge_status = rep('w', nrow(inc_mat))
  edge_status[membership_head != membership_tail] = 'b'
  return(edge_status)
}

# function to standardize Y
standardize <- function(x) {
  xmean = mean(x)
  x = x - xmean
  xscale = 2 * max(abs(x))
  x = x / xscale
  param = c('mean' = xmean, 'scale' = xscale)
  return(list(x = x, std_par = param))
}

# function to unstandardize Y
unstandardize <- function(x, std_par, nomean = F, s2 = F) {
  if(s2) {
    x = x * std_par['scale'] ^ 2
  } else {
    x = x * std_par['scale']
  }
  if(!nomean) x = x + std_par['mean']
  return(x)
}

# function to compute log-sum-exp
logSumExp <- function(x) {
  x_max = max(x)
  return(x_max + log( sum(exp(x - x_max)) ))
}


dentrigraph <- function(coords, threshold=1000){
  triangulation <- deldir(coords[,1], coords[,2])
  distances <- abs(triangulation$delsgs$x1 - triangulation$delsgs$x2) +
    abs(triangulation$delsgs$y1 - triangulation$delsgs$y2)
  #remove edges that are greater than thredhold
  edge.list <- as.matrix(triangulation$delsgs[distances < threshold, 5:6])
  
  if (length(unique(c(edge.list))) != nrow(coords))
    stop("Disconnected graph; try larger threshold")
  
  graph0 <- graph_from_edgelist(edge.list, directed = F)
  E(graph0)$weight = distances[distances < threshold]
  return(graph0)
}

# function to find distance between two spanning-tree partitions 
# i.e., number of edges with different status
partDist <- function(estatus_1, estatus_2) {
  return(sum(estatus_1 != estatus_2))
}

plotGraph <- function(coords, graph, title = NULL){
  edgelist = get.edgelist(graph, names = F) 
  edgedata = data.frame(coords[edgelist[, 1], ],coords[edgelist[, 2], ])
  colnames(edgedata) <- c("X1","Y1","X2","Y2")
  
  ggplot() + 
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2), data = edgedata, size = 0.5, colour = "grey")+
    labs(title = title, x = "lon", y = "lat")+
    theme(plot.title = element_text(hjust = 0.5))
}

plotGraphData <- function(coords, graph, Data, title = NULL, col_lim = NULL) {
  edgelist = get.edgelist(graph) 
  edgedata = data.frame(coords[edgelist[,1],], coords[edgelist[,2],])
  colnames(edgedata) <- c("X1","Y1","X2","Y2")
  if(missing(col_lim)) {col_lim = NA}
  
  ggplot() + 
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2), data = edgedata, size = 0.5, colour = "grey") + 
    geom_point(data = data.frame(coords), aes(lon, lat, color = Data))+
    scale_color_gradientn(colours = rainbow(10), limits = col_lim) +
    ggtitle(title) +
    theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
}
