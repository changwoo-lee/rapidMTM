#' Bayesian variable selection model based on Yang, Wainwright and Jordan, 2016, Annals of Statistics.

#' @param y n by 1, response vector.
#' @param X n by p, design matrix.
#' @param g positive real, hyperparameter 
#' @param kappa positive real, hyperparameter
#' @param s0 positive integer, hyperparameter (s_max in the manuscript)
#' @param burn integer, burn-in period
#' @param nmc integer, number of mcmc iteration after burn in
#' @param thin integer, thin-in
#' @param gammainit initial state gamma_0
#' @param verbose printing option
#' @param debug debugging option
#' @param preprocessed list of XtX, Xty and yty 
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
bvs_singletry <- function(y, X, g, kappa, s0, burn = 1000,
                          nmc = 5000, thin = 1, gammainit = NULL,
                          verbose = F, debug = F, preprocessed = NULL, 
                          truegamma=NULL){
  # for hitting the true
  found = F
  hit_iter = Inf
  hit_time = Inf
  truegammaidx = which(truegamma==1)
  
  p = ncol(X)
  n = nrow(X)
  gamma = numeric(p)
  if(is.null(gammainit)){
    gamma[sample(1:p, size = floor(s0/2))] <- 1 # random initialize gamma with size s0/2
  }else{
    gamma = gammainit
  }
  gammaidx = which(gamma == 1)
  N=burn+nmc
  effsamp=(N-burn)/thin
  ## output ##
  #betaout=matrix(0,p,effsamp)
  gammaout = list()
  accratioout = numeric(effsamp)
  logpostout = numeric(effsamp)
  Rsqout = numeric(effsamp)
  #sigmaSqout=rep(0,effsamp)
  #precal
  if(is.null(preprocessed)){
    cat("preprocessing XtX, Xty and yty...")
    Xtyoriginal = crossprod(X,y)
    #XtX = as(crossprod(X), "dspMatrix")
    XtX = crossprod(X)
    yty = sum(y^2)
    cat("done\n")
  }else{
    Xtyoriginal = preprocessed$Xty
    #XtX = as(preprocessed$XtX, "dspMatrix")
    XtX = preprocessed$XtX
    yty = as.numeric(preprocessed$yty)
  }
  if(length(gammaidx)>0){
    Xg = X[,gammaidx, drop = F]
    s = ncol(Xg)
    R = chol(crossprod(Xg))
  }else{
    s = 0
    R = matrix(numeric(0), 0, 0)
  }
  g.over.gp1 = g/(g+1)
  SSRout.current <- SSRgamma_fast(R, g.over.gp1, Xtyoriginal[gammaidx], yty)
  logSSR.current = log(SSRout.current$SSR)
  Rsq = SSRout.current$Rsq
  # MCMC time measure starts 
  mlikelihood_constant = lgamma(n/2) - n*log(pi)/2
  
  t_start = Sys.time()
  record.mcmc = 0
  for(imcmc in 1:N){
    #imcmc = imcmc + 1 
    acc = 0
    #if(runif(1)<0.5){# single flip, do nothing if gamma.prime is out of prior's support 
      idx = sample.int(p, size = 1) # this is scalable(not sensitive with p) 
      if(debug) cat(paste("single filp, s:",s,"\n"))
      # propose Rnew, snew, and gammaidxnew
      if(idx %in% gammaidx){# 1 if delete, 0 if add
        # delete variable
        snew = s - 1
        j = sum(gammaidx <= idx) # can be 1 to s
        Rnew = cholDel2(R, j)# delete one column
        gammaidxnew = gammaidx[-j]
      }else{ # add column
        snew = s + 1
        #if(snew == 1) browser()
        if(snew > s0){#do nothing if snew > s0
          record.mcmc = record.mcmc + 1
          gammaout[[record.mcmc]] <- gammaidx
          accratioout[record.mcmc] <- acc
          #log marginal likelihood
          logpost = mlikelihood_constant - (s/2)*log(1+g)  -s*kappa*log(p) - n/2*logSSR.current
          logpostout[record.mcmc] <- logpost
          next;
        }   
        j = sum(gammaidx <= idx) + 1 # CAN BE 1 TO s+1
        A2 = XtX[idx, idx, drop = F]
        A13 = XtX[gammaidx, idx, drop = F]
        if(j == 1) A1 = NULL else A1 = A13[1:(j-1),,drop = F]
        if(j == s+1) A3 = NULL else A3 = A13[j:s,,drop = F]
        Rnew = cholAdd2(R, j, A2, A1, A3)
        gammaidxnew = append(gammaidx, idx, after = j - 1)
      }
      Xtyselected <- Xtyoriginal[gammaidxnew]
      SSRout.proposal <- SSRgamma_fast(Rnew, g.over.gp1, Xtyselected, yty)
      logSSR.proposal = log(SSRout.proposal$SSR)
      logaccprob = (s - snew)*log(p^kappa*sqrt(1+g)) + (n/2)*(logSSR.current - logSSR.proposal)
      if(log(runif(1))<logaccprob){
        if(debug) cat(paste("single filp accepted, snew:",snew,"\n"))
        acc = 1
        #gamma[idx] <- 1 - gamma[idx]
        gammaidx = gammaidxnew
        s = snew
        R = Rnew
        logSSR.current = logSSR.proposal
        Rsq = SSRout.proposal$Rsq
      }
      
    #}else{ # double flip. do nothing when gamma is all zero 
    #   if(debug) cat(paste("double filp, s:",s,"\n"))
    #   if(s == 0){# null model, double filp not possible. Do nothing and go to next 
    #     record.mcmc = record.mcmc + 1
    #     gammaout[[record.mcmc]] <- gammaidx
    #     accratioout[record.mcmc] <- acc
    #     #log marginal likelihood
    #     logpost = mlikelihood_constant - (s/2)*log(1+g) -s*kappa*log(p) - n/2*logSSR.current
    #     logpostout[record.mcmc] <- logpost
    #     next;
    #   } 
    #   kidx = gammaidx[sample.int(s, size = 1)] # this is scalable(not sensitive with p) 
    #   #lidx.temp = sample.int(p-s, size = 1) # this is also scalable
    #   lidx.temp = sample.int(p-s, size = 1)
    #   #lidx = which(gamma==0)[lidx.temp] # this is NOT scalable; complexity linarly grows with p since length(gamma) = p
    #   #lidx.temp = (1:p)[-gammaidx]
    #   
    #   #alternatively... 
    #   lidx.tempnew = lidx.temp + sum(gammaidx <= lidx.temp)
    #   increment = sum(gammaidx <= lidx.tempnew & gammaidx > lidx.temp)
    #   while(increment > 0){ # while loop goes at most s0 times...
    #     lidx.temp = lidx.tempnew
    #     lidx.tempnew = lidx.temp + increment
    #     increment = sum(gammaidx <= lidx.tempnew & gammaidx > lidx.temp)
    #   }
    #   lidx = lidx.tempnew
    #   # 
    #   #lidx = lidx.temp + sum(gammaidx <= lidx.temp) # index l
    #   # delete kidx
    #   #browser()
    #   j = sum(gammaidx <= kidx)
    #   gammaidxnew = gammaidx[-j] 
    #   Rnew = cholDel2(R, j)
    #   snew = s-1
    #   
    #   # add lidx
    #   j = sum(gammaidxnew <= lidx) + 1 # CAN BE 1 TO s
    #   A2 = XtX[lidx, lidx, drop = F]
    #   A13 = XtX[gammaidxnew, lidx, drop = F]
    #   if(j == 1) A1 = NULL else A1 = A13[1:(j-1),,drop = F]
    #   if(j == snew+1) A3 = NULL else A3 = A13[j:snew,,drop = F]
    #   snew = snew+1 # 
    #   Rnew = cholAdd2(Rnew, j, A2, A1, A3)
    #   gammaidxnewnew = append(gammaidxnew, lidx, after = sum(gammaidxnew <= lidx)) # add lidx
    #   Xtyselected <- Xtyoriginal[gammaidxnewnew]
    #   SSRout.proposal <- SSRgamma_fast(Rnew, g.over.gp1, Xtyselected, yty)
    #   logSSR.proposal = log(SSRout.proposal$SSR)
    #   logaccprob = (n/2)*(logSSR.current - logSSR.proposal)
    #   
    #   #logSSR.proposal = log(SSRgamma(y, Xoriginal = X, g, gamma.prime))
    #   #logaccprob = (sum(gamma) - sum(gamma.prime))*log(p^kappa*sqrt(1+g)) + (n/2)*(logSSR.current - logSSR.proposal)
    #   if(log(runif(1))<logaccprob){
    #     if(debug) cat(paste("double filp accepted, snew:",snew,"\n"))
    #     acc = 1  
    #     # gamma[kidx] <- 1 - gamma[kidx]
    #     # gamma[lidx] <- 1 - gamma[lidx]
    #     gammaidx = gammaidxnewnew
    #     s = snew
    #     R = Rnew
    #     logSSR.current = logSSR.proposal
    #     Rsq = SSRout.proposal$Rsq
    #   }
    # }
    #if(identical(gammaidx,1:10)) return(imcmc)
    #step 1b(optional)draw phi(precision) and step 2(optional)draw beta: omitted
    if(verbose > 0) if(imcmc %% verbose == 0) {cat(paste("iteration",imcmc,"model size:",s,"\n"))};
    if(imcmc > burn && imcmc%%thin == 0)
    {
      record.mcmc = record.mcmc + 1
      gammaout[[record.mcmc]] <- gammaidx
      accratioout[record.mcmc] <- acc
      #log marginal likelihood
      logpost = mlikelihood_constant - (s/2)*log(1+g) -s*kappa*log(p) - n/2*logSSR.current
      logpostout[record.mcmc] <- logpost
      Rsqout[record.mcmc] <- Rsq
      
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
  cat("Elapsed time for",N,"MCMC iteration: ",mcmctime,"secs\n")
  
  bin <- tabulate(unlist(gammaout), nbins = p)
  PIP = bin/record.mcmc
  list(pip = PIP,
       gammaout = gammaout,
       accratio = mean(accratioout),
       logpostout = logpostout,
       Rsqout = Rsqout, 
       mcmctime = mcmctime,
       hit_time = hit_time,
       hit_iter = hit_iter)
}


library(mgcv)
library(Matrix)

# auxiliary functions
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


#helper function to compute the marginal probability
SSRgamma_fast <- function(R, g.over.gp1, Xty, yty){
  if(nrow(R)==0){
    Rsq = 0
    SSR = yty # null model
  }else{
    Rsq.numer = sum(forwardsolve(t(R),Xty)^2)
    SSR = yty - g.over.gp1 * Rsq.numer
    Rsq = Rsq.numer / yty
  }
  list(SSR = SSR, Rsq = Rsq)
}

