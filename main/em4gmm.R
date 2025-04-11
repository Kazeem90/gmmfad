# This code contains the GMMFAD algorithm and some associated functions. 
#
# It also contains the wrapper code for running the EMMIXmfa model under the same starting 
#  and stopping conditions as those of GMMFAD and the same initial parameters.
# It also contains a well optimized loglikelihood function for the MFA model.



library(MixSim)
library(dplyr)
library(stats)
library(fad)

#################################################################
###################################################################
#############        GMMFAD                #########################
###################################################################
#################################################################
## first, the pdf functions
logpdf<- function(X, mu, lambda, psi){
  n= dim(X)[1]
  p=length(psi)
  q=ifelse(is.vector(lambda), 1,  dim(lambda)[2] )
  psilambda <- 1/psi*lambda 
  M <- diag(q) + crossprod(lambda, psilambda) 
  invM <- chol2inv(chol(M))             
  invsig <- diag(1/psi)- tcrossprod( tcrossprod(psilambda, invM), psilambda) 
  logdetsig <-  sum(log(psi)) - log(det(invM)) 
  logconst  <- -0.5 * p * log(2*pi) - 0.5 * logdetsig
  X <- sweep(X, 2, mu)
  quadratic_term <-  rowSums(tcrossprod( X , invsig) * X)
  #####
  lpdf<-  logconst  - 0.5 * quadratic_term
  return(lpdf) 
}
###############################################################

############################
# gmm.fa2.k: this uses kmeans initialization. The most current version 
############################

library(fad)
gmm.fa.k <- function(X, K, q , tol = 1e-6, maxiter = 500, nstart= 50){
  X= as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  #kmeans
  kmeans_result <- kmeans(X, centers = K, nstart = 50)
  mu <- kmeans_result$centers
  sigma <- lapply(1:K, function(k) {
    cluster_data <- X[kmeans_result$cluster == k, ]
    cov(cluster_data)
  })
  sigma <- array(unlist(sigma), dim = c(p,p,K)) # to coerce sigma into an array.
  w <- as.vector(table(kmeans_result$cluster)) / n
  clusters <- kmeans_result$cluster
  n.clusters <- table(clusters)
  n.clusters <- as.vector(n.clusters)
  lambda <- array(dim = c(p, q, K) )
  psi <- matrix(nrow =K, ncol= p  )
  for (k in 1:K){ #cov2cor here might not be necessary.
    sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
    vars <- diag(sigma[,,k])
    start <- 1/ vars
    fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q, start= start, maxit = 500, tol = tol)) 
    sd = sqrt(vars)
    lambda[,,k] <-  lmd <-  fa$loadings 
    lambda[,,k]  <- sd * lmd
    psi[k,] <- ps <-  fa$uniquenesses
    psi[k,] <-  sd^2*ps
  }
  
  initials <- list(mean = mu, Sigma =sigma, weights = w , lambda =lambda, psi = psi)
  
  
  llh <- - Inf
  gamma <-  matrix(0, n, K)
  for (t in 1:maxiter){
    # E - Step
    # E - Step
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[,,k], psi[k,]))
    logp  <- sweep(logp, 2, log(w), "+")
    maxvals <- apply(logp, 1, max)
    logp  <- sweep(logp, 1, maxvals, "-")
    wprob <- exp(logp)
    sums  <- rowSums(wprob)
    gamma <- sweep(wprob, 1, sums , "/")
    
    # CM - Step 1
    gsum <- colSums(gamma)
    w <- gsum / n
    mu <- crossprod(gamma, X) / gsum
    for (k in 1:K){
      centered_X <- sweep(X, 2, mu[k,], "-")  # X - mu[k,]
      weighted_centered_X <- gamma[,k] * centered_X
      sigma[,,k] <- crossprod(weighted_centered_X, centered_X) / gsum[k]
      # CM - Step 2: compute lambdas and psis
      sigmaR <- suppressWarnings( cov2cor(sigma[,,k]) )
      vars <- diag(sigma[,,k])
      start <- 1/ vars
      fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q, start= start, maxit = 500, tol = tol)) 
      sd = sqrt(vars)
      lambda[,,k] <-  lam <-  fa$loadings 
      lambda[,,k] <- sd * lam # in covariance scale
      psi[k,] <- ps <-  fa$uniquenesses   # We can improve convergence speed by using the updated lambda to compute psi. 
      #First update sigma[,,k]= lambda[,,k]%*%lambda[,,k] + psi[k,] for each k. 
      psi[k,] <- sd^2*ps
    }
    llh0 <- llh
    #####################
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[,,k], psi[k,]))
    logp  <- sweep(logp, 2, log(w), "+")
    maxvals <- apply(logp, 1, max)
    logp  <- sweep(logp, 1, maxvals, "-")
    wprob <- exp(logp)
    sums  <- rowSums(wprob)
    llh   <- sum(maxvals, log(sums))
    # stopping criterion. 
    if  ( abs(llh-llh0) < tol * abs(llh) ){
      break
    }
  }
  gamma <- sweep(wprob, 1, sums , "/")
  clusters <- apply(gamma, 1, function(x) which.max(x) )
  iter <- t
  clusters <- apply(gamma, 1, function(x) which.max(x) )
  iter <- t
  #nparam <- K-1 + K * p + K * p* (p+1) / 2 + K*(p*q + p) # number of parameters 
  #nparam <- (K - 1) + K * p + p + K * (p * q - q * (q - 1)/2)
  nparam <-  (K - 1) + 2 * K * p + K * (p * q - q * (q - 1)/2)
  BIC <- -2*llh +  nparam * log(n) 
  
  analyser = list(clusters = clusters,  weights = w, 
                  lambda= lambda, psi=psi, sigma = sigma, niter =t,
                  means = mu, loglik= llh, BIC = BIC,  
                  nfactors = q, converged = t < maxiter, , gamma= gamma,
                  initial_param = initials
  )
  return(analyser)
}  
############################


############################
# gmm.fa: this uses random initialization.  
############################
library(fad)
gmm.fa <- function(X, K, q , tol =1e-6, maxiter = 500, nstart= 100, init=NULL){
  X= as.matrix(X)
  n= nrow(X)
  p= ncol(X)
  if (is.null(init)) {
    set.seed(12321)
  } else {
    set.seed(init)
  }
  
  mu <- matrix(rnorm(K*p),K )
  lambda <- array(rnorm(p*q*K), dim = c(p, q, K) )
  psi <- matrix(runif(K*p, 0.3,0.9) , nrow =K, ncol= p  )  
  sigma <- array(dim = c(p, p, K) )
  w<- rep(1/K, K)
  #for (k in 1:K){    
  #lambda[,,k] <- matrix(rnorm(p*q), p, q)
  #psi[k,] <- runif(p, 0.3,0.9)
  # sigma[,,k]<-  tcrossprod(lambda[,,k]) + diag(psi[k,]) 
  #}
  #pr<- abs(rnorm(K)) 
  #w<- rep(1/K, K) # pr/ sum(pr) 
  
  initials <- list(mean = mu, Sigma =sigma, weights = w , lambda =lambda, psi = psi)
  
  llh <- - Inf
  gamma <-  matrix(0, n, K)
  for (t in 1:maxiter){
    # E - Step
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[,,k], psi[k,]))
    logp  <- sweep(logp, 2, log(w), "+")
    maxvals <- apply(logp, 1, max)
    logp  <- sweep(logp, 1, maxvals, "-")
    wprob <- exp(logp)
    sums  <- rowSums(wprob)
    gamma <- sweep(wprob, 1, sums , "/")
    
    # CM - Step 1
    gsum <- colSums(gamma)
    w <- gsum / n
    mu <- crossprod(gamma, X) / gsum
    # CM - Step 2: compute lambdas and psis
    for (k in 1:K){
      centered_X <- sweep(X, 2, mu[k,], "-")  
      weighted_centered_X <- gamma[,k] * centered_X
      sigma[,,k] <- crossprod(weighted_centered_X, centered_X) / gsum[k]
      
      sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
      vars <- diag(sigma[,,k])
      start <- 1 / vars
      fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q, start= start, maxit = 500, tol = tol)) 
      sd = sqrt(vars)
      lambda[,,k] <-  lam <-  fa$loadings 
      lambda[,,k] <- sd * lam 
      psi[k,] <- ps <-  fa$uniquenesses   
      psi[k,] <- sd^2*ps
    }
    llh0 <- llh
    #####################
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[,,k], psi[k,]))
    logp  <- sweep(logp, 2, log(w), "+")
    maxvals <- apply(logp, 1, max)
    logp  <- sweep(logp, 1, maxvals, "-")
    wprob <- exp(logp)
    sums  <- rowSums(wprob)
    llh   <- sum(maxvals, log(sums))
    # stopping criterion. 
    if  ( abs(llh-llh0) < tol * abs(llh) ){
      break
    }
  }
  
  gamma <- sweep(wprob, 1, sums , "/")
  clusters <- apply(gamma, 1, function(x) which.max(x) )
  iter <- t
  
  nparam <-  (K - 1) + 2 * K * p + K * (p * q - q * (q - 1)/2)
  BIC <- -2*llh +  nparam * log(n) 
  
  analyser = list(clusters = clusters,  weights = w, 
                  lambda= lambda, psi=psi, sigma = sigma, niter =t,
                  means = mu, loglik= llh, BIC = BIC, 
                  nfactors = q, converged = t < maxiter, gamma= gamma, 
                  initial_param = initials
  )
  return(analyser)
}

##################################################################################



##################
# GMM wrapper: GMM with multiple initials 
##################
gmm.fad <- function(Y, K, q , tol = 1e-6, maxiter = 500, nstart= 20, 
                    init = NULL, innerNStart= 5){
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  llhs <- numeric( (nstart) )  
  res0 <- try({gmm.fa.k(Y, K, q, tol =tol,  maxiter = 1000) }, TRUE)
  llhs0 <- tryCatch({res0$loglik} , error= function(e) -Inf ) 
  
  if (is.null(init)) {
    set.seed(123456)
  } else {
    set.seed(init)
  }  
  seed <- runif(nstart)*10^8
  for (j in 1:nstart){
    res <- try({ gmm.fa(Y, K, q, tol =tol, init= seed[j], maxiter = 5) }, TRUE)
    llhs[j] <- tryCatch({res$loglik},  error= function(e) -Inf )     
  }
  llhs1 <- c(llhs, llhs0)
  if (sum(llhs1==-Inf)== (nstart + 1)) stop("All initial values and kmeans encountered errors. 
                      Number of factors or number of clusters might have 
                        been chosen too high compared to sample size
                        Also, try standardizing data. \n")
  irank <- order(llhs, decreasing = TRUE)[1:min(innerNStart, max(1,floor(nstart/4)))]
  
  for (j in irank){
    res2 <- try({
      gmm.fa(Y, K, q,  tol =tol, init= seed[j], maxiter = maxiter) 
    }, TRUE)
    llhs[j] <- tryCatch({res2$loglik},  error= function(e) -Inf)     
  }
  llhs1[1:nstart] <- llhs
  bestindex<- which.max(llhs1)
  if (bestindex== (nstart + 1) ){
    return( res0 )
  } else {
    res.f <-  gmm.fa(Y, K, q, tol = tol, init = seed[bestindex],  maxiter = maxiter)
    res.f$init = seed[bestindex]
    return( res.f )
  }
}

################################################################################
get_init_para<- function(p,g,q,init=NULL){
  K=g
  if (is.null(init)) {
    set.seed(12321)
  } else {
    set.seed(init)
  }
  
  mu <- matrix(rnorm(K*p),K )
  lambda <- array(rnorm(p*q*K), dim = c(p, q, K) )
  psi <- matrix(runif(K*p, 0.3,0.9) , nrow =K, ncol= p  )  
  w<- rep(1/K, K)
  psiarray <- array(dim=c(p,p, K))
  for (k in 1:K){psiarray[,,k]<- diag(psi[k,])}
  init_para <- list( g=K, q=q, pivec = w, mu  = t(mu), B = lambda, D=psiarray,  
                     sigma_type = "unique", D_type = "unique")
  return(init_para)
}

emmix_emEM <- function(Y, K, q, tol = 1e-6, maxiter = 500, nstart= 20, init = NULL, innerNStart= 5){
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  llhs <- numeric( (nstart) )  
  # llhs0 <- tryCatch({res0$loglik} , error= function(e) -Inf ) 
  
  if (is.null(init)) {
    set.seed(123456)
  } else {
    set.seed(init)
  }  
  seed <- runif(nstart)*10^8
  for (j in 1:nstart){
    init_para = get_init_para(p,K,q,init=seed[j])
    res <- mfa(Y, g=K, q=q, sigma_type = "unique", D_type = "unique", tol = tol, 
               conv_measure = 'ratio', itmax = 5, init_para= init_para)
    llhs[j] <- tryCatch({res$logL},  error= function(e) -Inf )     
  }
  #llhs1 <- c(llhs, llhs0)
  if (sum(llhs==-Inf)== (nstart)) stop("All initial values encountered errors. 
                      Number of factors or number of clusters might have 
                        been chosen too high compared to sample size
                        Also, try standardizing data. \n")
  irank <- order(llhs, decreasing = TRUE)[1:min(innerNStart, max(1,floor(nstart/4)))]
  for (j in irank){
    init_para = get_init_para(p,K,q,init=seed[j])
    res2 <- mfa(Y, g=K, q=q, sigma_type = "unique", D_type = "unique", tol = tol, 
                conv_measure = 'ratio', itmax = maxiter, init_para= init_para)
    llhs[j] <- tryCatch({res2$logL},  error= function(e) -Inf )      
  }
  bestindex<- which.max(llhs)
  init_para = get_init_para(p,K,q,init=seed[bestindex])
  res.f <- mfa(Y, g=K, q=q, sigma_type = "unique", D_type = "unique", tol = tol, 
               conv_measure = 'ratio', itmax = maxiter, init_para= init_para)
  #res.f$init = seed[bestindex]
  return( res.f )
}
#############################