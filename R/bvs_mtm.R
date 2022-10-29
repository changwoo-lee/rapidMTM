#' BVS with multiple-try metropolis with sqrt(u) weight function
#' Bayesian variable selection model based on Yang, Wainwright and Jordan (2016), On the computational complexity of high-dimensional Bayesian variable selection ,Annals of Statistics.
#'
#' @param y n by 1, response vector.
#' @param X n by p, design matrix.
#' @param g positive real, hyperparameter 
#' @param kappa positive real, hyperparameter
#' @param smax positive integer, hyperparameter (s_max in the manuscript)
#' @param burn integer, burn-in period
#' @param nmc integer, number of mcmc iteration after burn in
#' @param thin integer, thin-in
#' @param gammainit initial state gamma_0
#' @param verbose printing option
#' @param debug debugging option
#' @param preprocessed list of XtX, Xty and yty 
#' @param ntry the number of multiple-trials
#' @param truegamma (to measure the hit H) highest posterior prob. state
#'
#' @return list of:
#' pip = posterior inclusion prob.,
#' gammaout = list of the states (selected variables),
#' accratio = acceptance ratio,
#' logpostout = log-posterior,
#' Rsqout = R-squared, 
#' mcmctime = time took during mcmc,
#' hit_time = hit_time(T_H),
#' hit_iter = hit_iter(H).
#' @export
#'
#' @examples
bvs_mtm <- function(y, X, ntry, balancingft = NULL,
                         preprocessed = NULL,
                         burn = 1000, nmc = 5000, thin = 1, verbose = 1000,
                         gammainit = NULL, truegamma = NULL,
                         g = NULL, kappa =NULL, smax = NULL){
  ## balancing function
  if(!is.null(balancingft)){
    if(balancingft != "sqrt" & balancingft != "min" & balancingft != "max" ){
      stop("balancingft should be either `sqrt`, `min`, or `max`")
    }
  }else{
    print("balancingft is not provided, using ordinary weight function w(y|x) = p(y)")
    balancingft = "ord"
  }
  
  ### hyperparameter settings ###
  if(is.null(g)){
    g = p^3
  }
  #print(paste0("hyperparameter g is set as p^3 = ",g))
  g.over.gp1 = g/(g+1)
  
  if(is.null(kappa)){
    kappa = 2
  }
  #print(paste0("hyperparameter kappa is set as ",kappa))
  
  if(is.null(smax)){
    smax = 100
  }
  #print(paste0("hyperparameter smax is set as ",smax))
  
  y = as.vector(y)
  X = as.matrix(X)
  p = ncol(X)
  n = nrow(X)
  if(length(y)!=n) stop("length of y and number of rows of X does not match!")
  
  # precalculation
  if(is.null(preprocessed)){
    cat("preprocessing XtX, Xty and yty...")
    Xtyoriginal = crossprod(X,y)
    XtX = crossprod(X)
    yty = sum(y^2)
    cat("done\n")
  }else{
    Xtyoriginal = preprocessed$Xty
    XtX = preprocessed$XtX
    yty = as.numeric(preprocessed$yty)
  }
  diagXtX <- diag(XtX)
  
  # gamma
  if(is.null(gammainit)){
    gamma = numeric(p)
    gamma[sample(1:p, size = floor(smax/2))] <- 1 # random initialize gamma with size smax/2
  }else{
    gamma = gammainit
  }
  gammaidx = which(gamma==1)
  
  # true gamma, if any
  found = F
  hit_iter = Inf
  hit_time = Inf
  if(!is.null(truegamma)){
    if(length(truegamma)!=p) stop("truegamma length and number of columns of X does not match")
    truegammaidx = which(truegamma==1)
    cat(paste0("true gamma provided with model size ",length(truegammaidx),", measuring hitting iteration and hitting (wall-clock) time...\n"))
  } 
  
  # MCMC 
  nmcmc = burn+nmc
  nsamples = floor((nmcmc-burn)/thin)
  
  
  ## output ##
  gammaout = list()
  accratioout = numeric(nsamples)
  logpostout = numeric(nsamples)
  
  # s: size of gamma
  # R: chol(X_gamma^T X_gamma) such that R^T R = X_gamma^T X_gamma
  if(length(gammaidx)>0){
    Xg = X[,gammaidx, drop = F]
    s = ncol(Xg)
    R = chol(crossprod(Xg)) 
  }else{
    s = 0
    R = matrix(numeric(0), 0, 0)
  }
  # bsol = (R^T)^(-1)* X_gamma^T * y so that bsol^T * bsol = y^T X_gamma (X_gamma^T X_gamma)^(-1) X_gamma^T y
  if(nrow(R)==0){
    SSR = yty # null model
  }else{
    bsol = backsolve(R, Xtyoriginal[gammaidx], transpose = T)
    SSR = yty - g.over.gp1 * sum(bsol^2)
  }
  
  # log posterior ( see eq.14, except constant C )
  SSR.current = SSR
  likelihood_constant =  kappa*log(p) + log(1+g)/2 
  logpost.current = -s*likelihood_constant - (n/2)*log(SSR.current)
  
  ###  MCMC time measure starts  ####
  cat(paste("Run MCMC, initial model size:",s,"\n"))
  t_start = Sys.time()
  record.mcmc = 0
  for(imcmc in 1:nmcmc){
    #imcmc = imcmc + 1 
    acc = 0
    
    # if s < smax,
    
    # Step 1. ntry independent single flip proposals
    # sample the indicies that will change, not length p vector
    idx_K = sample.int(p, size = ntry, replace = T) # this is scalable (not sensitive with p)
    
    # Delete variables idx
    Kidx_drop = which(idx_K%in%gammaidx) # indicies of idx_K that is included in the gammaidx 
    # Add variables idx
    Kidx_add = which(!(idx_K%in%gammaidx))
    
    logweights = numeric(ntry)
    logweights_star = numeric(ntry-1)
    
    logpost_multiple = numeric(ntry) # placeholder to save log posterior probabilities for each ntry
    
    Rnew_list <- list()
    
    # calculate log weights for delete moves
    for(jj in Kidx_drop){
      j = sum(gammaidx <= idx_K[jj]) # can be 1 to s
      Rnew_list[[jj]] <- Rnew <- cholDel2(R, j)# delete one column
      gammaidxnew = gammaidx[-j]
      if(nrow(Rnew)==0){
        SSR = yty # null model
      }else{
        SSR = yty - g.over.gp1 * sum(backsolve(Rnew, Xtyoriginal[gammaidxnew], transpose = T)^2)
      }
      logpost_multiple[jj] = -(s-1)*likelihood_constant - (n/2)*log(SSR)
    }
    
    # calculate log weights for add moves, with matrix operation 
    if( s+1 > smax){#do nothing if snew > smax
      logpost_multiple[Kidx_add] = -Inf
    }else if(length(Kidx_add)==0){ # no add proposal 
      # do nothing
    }else{
      if(s > 0){
        jvec = idx_K[Kidx_add]
        n.kadd = length(Kidx_add)
        S12.parallel = matrix(backsolve(R, XtX[gammaidx, jvec, drop = F], transpose = T), ncol = n.kadd) # s x n.kadd matrix
        S22.parallel = sqrt(as.numeric(diagXtX[jvec] - colSums(S12.parallel^2))) # length n.kadd vector
        #if(any(is.na(S22.parallel))) browser()
        bsol_norm_diff = ((Xtyoriginal[jvec] - crossprod(S12.parallel,bsol))/S22.parallel)^2 # length n.kadd vector
        SSR.proposal = SSR.current - g.over.gp1*bsol_norm_diff # amount SSR decreasing, length n.kadd vector
        logpost_multiple[Kidx_add] = - (s+1)*likelihood_constant - (n/2)*log(SSR.proposal)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
        #browser()
      }else{ # add from null model
        jvec = idx_K[Kidx_add]
        # no S12
        S22.parallel = sqrt(as.numeric(diagXtX[jvec])) # length kadd vector
        bsol_norm_diff = (Xtyoriginal[jvec]/S22.parallel)^2
        SSR.proposal = SSR.current - g.over.gp1*bsol_norm_diff 
        logpost_multiple[Kidx_add] = - (s+1)*likelihood_constant - (n/2)*log(SSR.proposal)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
      }
    }
    
    if(balancingft == "sqrt"){
      logweights = (logpost_multiple - logpost.current)/2
    }else if(balancingft == "min"){
      logweights = pmin(logpost_multiple - logpost.current, 0)
    }else if(balancingft == "max"){
      logweights = pmax(logpost_multiple - logpost.current, 0)
    }else{
      logweights = logpost_multiple
    }
    accratio.numer = matrixStats::logSumExp(logweights) #log scale
    weights_normalized = exp(logweights -  accratio.numer)
    
    # Step 2. Select j with probability proportional to w
    proposed.Kidx = sample.int(ntry, size = 1, prob = weights_normalized) # final proposal index among 1:ntry
    proposed.idx = idx_K[proposed.Kidx] #FINAL PROPOSAL index
    
    
    #############################################
    if(proposed.Kidx %in% Kidx_drop){# drop is proposed
      
      j = sum(gammaidx <= proposed.idx)
      Rnew = Rnew_list[[proposed.Kidx]]
      gammaidxnew = gammaidx[-j]
      if(length(gammaidxnew)>0){
        bsolnew = backsolve(Rnew,Xtyoriginal[gammaidxnew], transpose = T)
        SSR.proposed = yty - g.over.gp1 * sum(bsolnew^2)
      }else{ # null model
        bsolnew = 0
        SSR.proposed = yty
      }
      snew = s-1
      
    }else{ # add is proposed
      
      j = sum(gammaidx <= proposed.idx) + 1 # CAN BE 1 TO s+1
      A2 = diagXtX[proposed.idx]
      A13 = XtX[gammaidx, proposed.idx, drop = F]
      if(j == 1) A1 = NULL else A1 = A13[1:(j-1),,drop = F]
      if(j == s+1) A3 = NULL else A3 = A13[j:s,,drop = F]
      Rnew = cholAdd2(R, j, A2, A1, A3)
      gammaidxnew = append(gammaidx, proposed.idx, after = j - 1)
      bsolnew = backsolve(Rnew,Xtyoriginal[gammaidxnew], transpose = T)
      SSR.proposed = yty - g.over.gp1 * sum(bsolnew^2)
      snew = s+1
    }
    
    logpost.proposed = logpost_multiple[proposed.Kidx]
    
    ### step 3 ##########################################
    
    idx_K_star = sample.int(p, size = ntry-1, replace = T) # auxiliary 
    
    Kidx_drop_star = which(idx_K_star%in%gammaidxnew) # indicies of idx_K that is included in the gammaidx 
    Kidx_add_star = which(!(idx_K_star%in%gammaidxnew))
    logpost_multiple_star = numeric(ntry-1) # length K -1, save log posterior probabilities
    
    
    # delete move
    
    for(jj in Kidx_drop_star){
      j = sum(gammaidxnew <= idx_K_star[jj]) # can be 1 to s
      Rnew_star <- cholDel2(Rnew, j)# delete one column
      gammaidxnew_star = gammaidxnew[-j]
      if(nrow(Rnew_star)==0){
        SSR_star = yty # null model
      }else{
        SSR_star = yty - g.over.gp1 * sum(backsolve(Rnew_star, Xtyoriginal[gammaidxnew_star], transpose = T)^2)
      }
      logpost_multiple_star[jj] =  - (snew-1)*likelihood_constant - (n/2)*log(SSR_star)
    }
    
    # add move, using matrix operation
    
    if( snew+1 > smax){#do nothing if snew > smax
      logpost_multiple_star[Kidx_add_star] = -Inf
      #browser()
      #K_weight_normalized[Kidx_add_star] = 0
    }else if(length(Kidx_add_star)==0){ # no add proposal 
      # do nothing
    }else{
      if(snew > 0){
        jvec = idx_K_star[Kidx_add_star]
        n.kadd = length(Kidx_add_star)
        S12.parallel = matrix(backsolve(Rnew, XtX[gammaidxnew, jvec, drop = F], transpose = T), ncol = n.kadd) # s x n.kadd matrix
        S22.parallel = sqrt(as.numeric(diagXtX[jvec] - colSums(S12.parallel^2))) # length n.kadd vector
        if(any(is.na(S22.parallel))) browser()
        bsol_norm_diff = ((Xtyoriginal[jvec] - crossprod(S12.parallel,bsolnew))/S22.parallel)^2 # length n.kadd vector
        
        SSR.proposal_star = SSR.proposed - g.over.gp1*bsol_norm_diff # amount SSR decreasing, length n.kadd vector
        logpost_multiple_star[Kidx_add_star] = - (snew+1)*likelihood_constant - (n/2)*log(SSR.proposal_star)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
        #browser()
      }else{ # add from null model
        jvec = idx_K_star[Kidx_add_star]
        # no S12
        S22.parallel = sqrt(as.numeric(diagXtX[jvec])) # length kadd vector
        bsol_norm_diff = (Xtyoriginal[jvec]/S22.parallel)^2
        SSR.proposal_star = SSR.proposed - g.over.gp1*bsol_norm_diff 
        logpost_multiple_star[Kidx_add_star] = - (snew+1)*likelihood_constant - (n/2)*log(SSR.proposal_star)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
      }
    }
    
    # w(x_t | y) and w(x_l^* | y)

    if(balancingft == "sqrt"){
      logweights_star = (c(logpost.current, logpost_multiple_star) - logpost.proposed)/2
    }else if(balancingft == "min"){
      logweights_star = pmin(c(logpost.current, logpost_multiple_star) - logpost.proposed, 0)
    }else if(balancingft == "max"){
      logweights_star = pmax(c(logpost.current, logpost_multiple_star) - logpost.proposed, 0)
    }else{
      logweights_star = c(logpost.current, logpost_multiple_star)
    }
    
    #########################################################################################
    #if(any(is.na(logpost_multiple_star))) browser()
    
    accratio.denom = matrixStats::logSumExp(logweights_star) # log scale
    
    if(log(runif(1)) < accratio.numer - accratio.denom){
      acc = 1
      #if(debug) cat(paste("single filp accepted,\n"))
      R = Rnew
      gammaidx = gammaidxnew
      bsol = bsolnew
      SSR.current = SSR.proposed
      s = snew
      logpost.current = logpost.proposed
    }
    
    
    #if(length(gammaidx)==9) return(imcmc)
    #step 1b(optional)draw phi(precision) and step 2(optional)draw beta: omitted
    if(imcmc %% verbose == 0) {cat(paste("iteration",imcmc,"model size:",s,"\n"))};
    if(imcmc > burn && imcmc%%thin == 0)
    {
      record.mcmc = record.mcmc + 1
      gammaout[[record.mcmc]] <- gammaidx
      accratioout[record.mcmc] <- acc
      logpostout[record.mcmc] <- logpost.current
      
      # record hit-time
      if(!is.null(truegamma) & !found){
        if(length(gammaidx)==length(truegammaidx)){
          if(all(sort(gammaidx) == truegammaidx)){
            hit_iter = imcmc
            found = T
            # summarize
            hit_time = difftime(Sys.time(), t_start, units = "secs")
          }
        }
      }
      
    }
  }
  # summarize
  mcmctime = difftime(Sys.time(), t_start, units = "secs")
  cat("Elapsed time for",nmcmc,"MCMC iteration: ",mcmctime,"secs\n")
  #browser()
  bin <- tabulate(unlist(gammaout), nbins = p)
  PIP = bin/record.mcmc
  
  list(pip = PIP,
       gammaout = gammaout,
       accratio = mean(accratioout),
       logpostout = logpostout,
       mcmctime = mcmctime,
       hit_time = hit_time,
       hit_iter = hit_iter)
}




library(mgcv)
library(Matrix)

cholDel2 <- function(R, j){
  n = nrow(R)
  if(j == n){
    R_new = R[-n, -n, drop = F]
  }else{
    R_new = mgcv::choldrop(R, j)
  }
  R_new
}



cholAdd2 <- function(R, j, A2, A1 = NULL, A3 = NULL) {
  n = nrow(R)
  if(j == n+1) {
    R_new = matrix(0, nrow = n+1, ncol = n+1)
    if(n > 0){
      R_new[1:n,1:n] = R
      S12 = drop(backsolve(R, A1, transpose = T))
      S22 = sqrt(as.numeric(A2 - sum(S12^2)))
      R_new[1:(n+1), n+1] = c(S12, S22)
    }else{
      R_new[1,1] = sqrt(as.numeric(A2)) ### fixed - 20210331 / 20210412
    }
  } else {
    R_new = matrix(0, nrow = n+1, ncol = n+1)
    if(j > 1) {
      R11 = R[1:(j-1), 1:(j-1), drop = F]
      R_new[1:(j-1), 1:(j-1)] = R11
      S12 = backsolve(R11, A1, transpose = T)
      R_new[1:(j-1), j] = S12
      S13 = R[1:(j-1), j:n, drop = F]
      R_new[1:(j-1), (j+1):(n+1)] = S13
      
      S22 = sqrt(as.numeric(A2 - sum(S12^2)))
      R_new[j, j] = S22
      S23 = (t(A3) - crossprod(S12, S13)) / S22
      R_new[j, (j+1):(n+1)] = S23
      S33 = mgcv::cholup(R[j:n, j:n, drop = F], S23, up = FALSE) # downdate
      R_new[(j+1):(n+1), (j+1):(n+1)] = S33
    } else {
      S22 = sqrt(as.numeric(A2))
      R_new[1, 1] = S22
      S23 = as.numeric(A3) / S22
      R_new[1, 2:(n+1)] = S23
      S33 = mgcv::cholup(R, S23, up = FALSE) # downdate
      R_new[2:(n+1), 2:(n+1)] = S33
    }
  }
  return(R_new)
}

