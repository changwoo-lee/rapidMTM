#' SBM with multiple-try metropolis with sqrt(u) weight function
#' Bayesian stochastic block model based on Zhuo and Gao (2021), Mixing Time of Metropolis-Hastings for Bayesian Community Detection, JMLR.
#'
#' @param A n by n (for notation in paper, p by p) binary adjacency matrix
#' @param k number of blocks, fixed and known
#' @param ntry number of trials for MTM
#' @param z_init initial state(partition), represented as a length n integer vector
#' @param N_iter number of mcmc iteration
#' @param a first parameter of beta prior, default 1 
#' @param b second parameter of beta prior, default 1
#' @param savez logical, whether save partitions z?
#' @param true_z (to measure the hit H) highest posterior prob. state
#' @param sparse logical, using sparse matrix multiplication?
#'
#' @return list of:
#' (if savez = T) z_post: partitions
#' loglik_post: log-posterior,
#' accept_post: acceptance ratio,
#' mcmctime = time took during mcmc,
#' hit_time = hit_time(T_H),
#' hit_iter = hit_iter(H).
#' @export
#'
#' @examples
sbm_mtm_sqrt<- function(A, k =5, ntry = 10, z_init = NULL, N_iter=10000, a=1, b=1, 
                                 savez = F, true_z = NULL, sparse = F){
    # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  # for hitting the true
  found = F
  hit_iter = Inf
  hit_time = Inf
  if(!is.null(true_z)){
    true_z = cluster_sort(true_z)
    true_nclust = tabulate(true_z)
    true_nclust_sorted = sort(true_nclust)
  }
  
  
  n <- nrow(A)
  if(!sparse){
    #mode(A) = "logical" # fast indexing
  }else{
    A = Matrix::Matrix(A)
    A = as(A, "nMatrix") ## important, improves speed almost x2
  }
  
  Adj_Edgelist = list()
  for(i in 1:n) Adj_Edgelist[[i]] = which(A[,i]==1)
    
  if(is.null(z_init)){
    z_init = rep(1:k, each = n/k, length = n)
    #z_init = sample(z_init)
    #z_init = rep(1:k, times = c(5, 16, 19, 26, 34))
  }else{
    if(length(z_init)!=n) stop("wrong z_init length")
    if(length(unique(z_init))!=k) stop("wrong ncluster of z_init") 
  }
  #z_init = GPPM::cluster_sort(z_init)
  z = z_init
  mode(z) = "integer"
  # cluster assignments are encoded in two equivalent waAs:
  # [i] a VxH matrix Z, s.t. Z[i,h]=1{node i in cluster h}, faster to use within each iteration
  Z <- onehot(z, n, k)

  # [ii] a vector of length n containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[i,t]=h if node i is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  if(savez) z_post <- matrix(NA,N_iter, n)
  loglik_post <- numeric(N_iter)
  accept_post <- numeric(N_iter)
  # Create the matrix with block connections
  temp   <- A%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # number of edges between clusters, represented by vector
  m_vector = m_full[upper.tri(m_full, diag = T)]
  
  nclust = tabulate(z)
  
  N_full = tcrossprod(nclust)
  diag(N_full) = nclust*(nclust-1)/2
  # number of possible edges between clusters, represented by vector
  N_vector = N_full[upper.tri(N_full, diag = T)]
  
  log_lik = sum(lbeta(m_vector + a, N_vector - m_vector + b))
  
  
  # Index lookuptable
  index_lookup = matrix(0L, k, k)
  index_lookup[upper.tri(index_lookup, diag = T)] <- 1:(k*(k+1)/2)
  index_lookup = (index_lookup + t(index_lookup))
  diag(index_lookup) = diag(index_lookup)/2
  #index_lookup
  appending_idx1 = rep((0:(ntry-1))*(k*(k+1)/2), each = k)
  appending_idx2 = (0:(ntry-1))*(k*(k+1)/2)
  
  appending_idx3 = rep((0:(ntry-2))*(k*(k+1)/2), each = k)
  appending_idx4 = (0:(ntry-2))*(k*(k+1)/2)
  
  # ----------------------------------------------
  # Beginning of the M-H sampler
  # ----------------------------------------------
  one_to_k = 1:k
  
  t_start = Sys.time()
  for (imcmc in 1:N_iter){
    
    # j is now vector 
    j = sample.int(n, size = ntry, replace = T)
    j_oldclust = z[j]
    # cluster number change is not allowed
    if(any(nclust[j_oldclust]==1)) next;
    
    # 2-multi. 
    rand_km1 = sample.int(k-1, size = ntry, replace = T)
    j_newclust = rand_km1 + (j_oldclust<=rand_km1)
    
    # 3. calculate likelihood ratio
    #browser()
    if(sparse){
      r_v = as.vector(crossprod(Z, A[,j])) # same as tabulate(z[which(A[,j])], nbins = k)
    }else{
      r_v = as.vector(vapply(Adj_Edgelist[j], FUN = myf, z = z, k = k, FUN.VALUE = integer(k)))
    }
    # microbenchmark(
    # #r_v = as.vector(crossprod(Z, A[,j])), # same as tabulate(z[which(A[,j])], nbins = k)
    # #r_v2 = as.vector(sapply(Adj_Edgelist[j], FUN = f, z = z)),
    # r_v3 = as.vector(vapply(Adj_Edgelist[j], FUN = myf, z = z, k=k,FUN.VALUE = integer(k))),
    # r_v4 = as.vector(vapply(Adj_Edgelist[j], FUN = myf, z = zint, k=k,FUN.VALUE = integer(k)))
    # )
    
    oldidx = as.vector(index_lookup[,j_oldclust]) + appending_idx1 # appended, vector
    newidx = as.vector(index_lookup[,j_newclust]) + appending_idx1 # appended, vector
    idx_plusone = index_lookup[cbind(j_oldclust,j_oldclust)] + appending_idx2
    idx_minusone = index_lookup[cbind(j_oldclust,j_newclust)]+ appending_idx2
    #browser()
    # each column is each try
    m_vector_appended_new = rep(m_vector, ntry)
    m_vector_appended_new[oldidx] = m_vector_appended_new[oldidx] - r_v
    m_vector_appended_new[newidx] = m_vector_appended_new[newidx] + r_v
    
    m_matrix_new = matrix(m_vector_appended_new, ncol = ntry)
    
    # N_vector
    N_vector_appended_new = rep(N_vector, ntry)
    N_vector_appended_new[oldidx] = N_vector_appended_new[oldidx] - nclust
    N_vector_appended_new[newidx] = N_vector_appended_new[newidx] + nclust
    N_vector_appended_new[idx_plusone] = N_vector_appended_new[idx_plusone] +1
    N_vector_appended_new[idx_minusone] = N_vector_appended_new[idx_minusone] -1
    
    N_matrix_new = matrix(N_vector_appended_new, ncol = ntry)
    
    # ntry vector
    log_lik_new = colSums(lbeta(m_matrix_new + a, N_matrix_new - m_matrix_new + b))
    
    # square root weighting
    
    logweights = (log_lik_new - log_lik)/2
    accratio.numer = matrixStats::logSumExp(logweights) #log scale
    weights_normalized = exp(logweights - accratio.numer)
    
    proposed.Kidx = sample.int(ntry, size = 1, prob = weights_normalized) # final proposal index of 1:K
    proposed.idx = j[proposed.Kidx] # proposed single flip
    
    # next 
    #browser()
    proposed.z = z
    proposed.z[proposed.idx] = j_newclust[proposed.Kidx]
    proposed.m_vector = m_matrix_new[,proposed.Kidx]
    proposed.N_vector = N_matrix_new[,proposed.Kidx]
    proposed.loglik = log_lik_new[proposed.Kidx]
    
    proposed.nclust = nclust
    proposed.nclust[j_newclust[proposed.Kidx]] = proposed.nclust[j_newclust[proposed.Kidx]] + 1
    proposed.nclust[j_oldclust[proposed.Kidx]] = proposed.nclust[j_oldclust[proposed.Kidx]] - 1
    
    if(sparse){
      proposed.Z = Z
      proposed.Z[j[proposed.Kidx], j_oldclust[proposed.Kidx]] = 0
      proposed.Z[j[proposed.Kidx], j_newclust[proposed.Kidx]] = 1
    }
    
    ## single flip from the proposed state
    jstar = sample.int(n, size = ntry -1 , replace = T)
    jstar_oldclust = proposed.z[jstar]
    
    # cluster number change is not allowed
    if(any(proposed.nclust[jstar_oldclust]==1)) next;
    
    # 2-multi. 
    rand_km1 = sample.int(k-1, size = ntry-1, replace = T)
    jstar_newclust = rand_km1 + (jstar_oldclust<=rand_km1)
    
    
    # 3. calculate likelihood ratio
    if(sparse){
      r_v = as.vector(crossprod(proposed.Z, A[,jstar]))
    }else{
      r_v = as.vector(vapply(Adj_Edgelist[jstar], FUN = myf, z = proposed.z, k = k, FUN.VALUE = integer(k)))
    }
    
    oldidx = as.vector(index_lookup[,jstar_oldclust]) + appending_idx3 # appended, vector
    newidx = as.vector(index_lookup[,jstar_newclust]) + appending_idx3 # appended, vector
    idx_plusone = index_lookup[cbind(jstar_oldclust,jstar_oldclust)] + appending_idx4
    idx_minusone = index_lookup[cbind(jstar_oldclust,jstar_newclust)]+ appending_idx4
    
    # each column is each try
    m_vector_appended_new = rep(proposed.m_vector, ntry - 1)
    m_vector_appended_new[oldidx] = m_vector_appended_new[oldidx] - r_v
    m_vector_appended_new[newidx] = m_vector_appended_new[newidx] + r_v
    
    m_matrix_newstar = matrix(m_vector_appended_new, ncol = ntry - 1)
    
    # N_vector
    N_vector_appended_new = rep(proposed.N_vector, ntry - 1)
    N_vector_appended_new[oldidx] = N_vector_appended_new[oldidx] - proposed.nclust
    N_vector_appended_new[newidx] = N_vector_appended_new[newidx] + proposed.nclust
    N_vector_appended_new[idx_plusone] = N_vector_appended_new[idx_plusone] +1
    N_vector_appended_new[idx_minusone] = N_vector_appended_new[idx_minusone] -1
    
    N_matrix_newstar = matrix(N_vector_appended_new, ncol = ntry - 1)
    
    # ntry -1 vector
    log_lik_newstar = colSums(lbeta(m_matrix_newstar + a, N_matrix_newstar - m_matrix_newstar + b))
    
    logweights_star = (log_lik_newstar - proposed.loglik)/2
    
    accratio.denom = matrixStats::logSumExp(c(logweights_star, -logweights[proposed.Kidx] )) # log scale
    
    accept = F
    if(log(runif(1)) < accratio.numer - accratio.denom){
      accept = T
      z = proposed.z
      if(sparse) Z = proposed.Z
      m_vector = proposed.m_vector 
      N_vector = proposed.N_vector
      log_lik = proposed.loglik
      nclust = proposed.nclust 
      
    }
    
    
    if(savez) z_post[imcmc,] <- z
    loglik_post[imcmc] <- log_lik
    accept_post[imcmc] <- accept
    
    # record hit-time
    if(!is.null(true_z) & !found){
      if(all(sort(nclust)==true_nclust_sorted)){ 
        if(all(cluster_sort(z) == true_z)){
          hit_iter = imcmc
          found = T
          # summarize
          hit_time = difftime(Sys.time(), t_start, units = "secs")
        }
      }
    }
    
    #print(table(z_post[,t])) 
    #if (imcmc%%1000 == 0){print(paste("Iteration:", imcmc))}
  }
  mcmctime = difftime(Sys.time(), t_start, units = "secs")
  cat("Elapsed time for",N_iter,"MCMC iteration: ",mcmctime,"secs\n")
  
  out = list()
  if(savez) out$z_post = z_post
  out$loglik_post = loglik_post
  out$accept_post = accept_post
  out$mcmctime = mcmctime
  out$hit_time = hit_time
  out$hit_iter = hit_iter
  return(out)
}





onehot <- function(z, n, k){
  Z <- matrix(0,n,k)
  for (j in 1:k){ # not i in 1:n
    Z[which(z==j),j] <- 1
  }
  return(Z)
}

cluster_sort <- function(x){ 
  nunique = length(unique(x))
  temp = factor(x, labels = 1:nunique)
  temp = as.numeric(as.character(temp))
  res = replace(temp, unique(temp), 1:nunique)[temp]
  res
}

calculatelog_lik <- function(A, z){
  n = length(z)
  k = length(unique(z))
  Z = onehot(z, n, k)
  # Create the matrix with block connections
  mode(A) = "logical"
  temp   <- A%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),k)
  
  # number of edges between clusters, represented by vector
  m_vector = m_full[upper.tri(m_full, diag = T)]
  
  nclust = tabulate(z)
  
  N_full = tcrossprod(nclust)
  diag(N_full) = nclust*(nclust-1)/2
  # number of possible edges between clusters, represented by vector
  N_vector = N_full[upper.tri(N_full, diag = T)]
  
  log_lik = sum(lbeta(m_vector + 1, N_vector - m_vector + 1))
  log_lik
}

require(Matrix)
myf = function(x, z, k){
  tabulate(z[x], nbins = k)
}


