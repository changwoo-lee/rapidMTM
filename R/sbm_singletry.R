#' Bayesian stochastic block model based on Zhuo and Gao (2021), Mixing Time of Metropolis-Hastings for Bayesian Community Detection, JMLR.
#'
#' @param A n by n (for notation in paper, p by p) binary adjacency matrix
#' @param k number of blocks, fixed and known
#' @param z_init initial state(partition), represented as a length n integer vector
#' @param N_iter number of mcmc iteration
#' @param a first parameter of beta prior, default 1 
#' @param b second parameter of beta prior, default 1
#' @param savez logical, whether save partitions z?
#' @param true_z (to measure the hit H) highest posterior prob. state
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
sbm_singletry <- function(A, k =5, z_init = NULL, N_iter=10000, a=1, b=1, savez = F, true_z = NULL){
  
  found = F
  hit_iter = Inf
  hit_time = Inf
  if(!is.null(true_z)){
    true_z = cluster_sort(true_z)
    true_nclust = tabulate(true_z)
    true_nclust_sorted = sort(true_nclust)
  }
  
  # singletry
  #ntry= 1
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  n <- nrow(A)
  mode(A) = "logical" # fast indexing

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
  #appending_idx1 = rep((0:(ntry-1))*(k*(k+1)/2), each = k)
  #appending_idx2 = (0:(ntry-1))*(k*(k+1)/2)

  # ----------------------------------------------
  # Beginning of the M-H sampler
  # ----------------------------------------------
  one_to_k = 1:k
  
  
  t_start = Sys.time()
  for (imcmc in 1:N_iter){
    # ntry = 1
    
    # 1. choose an index j uniformly at random
    j = sample.int(n, size = 1)
    j_oldclust = z[j]
    #browser()
    # cluster number change is not allowed
    if(nclust[j_oldclust]==1){
      print("singleton!")
      next; 
    }
    
    # 2. randomly assign new label to get new assignment
    j_newclust = sample(one_to_k[-j_oldclust], size =1)
    
    
    
    # 3. calculate likelihood ratio
    
    #r_v = crossprod(Z, A[,j]) # same as 
    r_v = tabulate(z[which(A[,j])], nbins = k)
    
    oldidx = index_lookup[,j_oldclust]
    newidx = index_lookup[,j_newclust]
    
    # m_vector
    m_vector_new = m_vector
    m_vector_new[oldidx] = m_vector_new[oldidx] - r_v
    m_vector_new[newidx] = m_vector_new[newidx] + r_v
    
    # N_vector
    N_vector_new = N_vector
    N_vector_new[oldidx] = N_vector_new[oldidx] - nclust
    N_vector_new[newidx] = N_vector_new[newidx] + nclust
    N_vector_new[index_lookup[j_oldclust,j_oldclust]] = N_vector_new[index_lookup[j_oldclust,j_oldclust]] + 1
    N_vector_new[index_lookup[j_oldclust,j_newclust]] = N_vector_new[index_lookup[j_oldclust,j_newclust]] - 1
    
    # comparison..
    # nclust_new = nclust
    # nclust_new[j_oldclust] = nclust_new[j_oldclust]-1
    # nclust_new[j_newclust] = nclust_new[j_newclust]+1
    # 
    # N_full_new = tcrossprod(nclust_new)
    # diag(N_full_new) = nclust_new*(nclust_new-1)/2
    # N_full_new
    
    log_lik_new = sum(lbeta(m_vector_new + a, N_vector_new - m_vector_new + b))
    
    accept = F
    if(log(runif(1)) < log_lik_new - log_lik){
      accept = T
      log_lik = log_lik_new  
      # update z
      
      z[j] <- j_newclust
      #z = GPPM::cluster_sort(z) don't need this
      # update Z (one-hot encoding matrix)
      #Z = onehot(z, n, k)
      
      # update nclust
      #nclust = tabulate(z)}, #maybe simplified
      nclust[j_newclust] = nclust[j_newclust] + 1
      nclust[j_oldclust] = nclust[j_oldclust] - 1
      
      # update m vector, n vector
      m_vector = m_vector_new
      N_vector = N_vector_new
    }
    #browser()
    # store cluster assignments at time imcmc in matrix z_post s.t.
    # z_post[i,t]=h if node i is in cluster h at iteration t
    if(savez) z_post[imcmc,] <- z
    loglik_post[imcmc] <- log_lik
    accept_post[imcmc] <- accept
    #print(table(z_post[,t])) 
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






