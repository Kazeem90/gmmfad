library(MixSim)
library(dplyr)
library(stats)
################################################################
## first, the pdf functions
logpdf<- function(X, mu, lambda, psi){
  n= dim(X)[1]
  p=length(psi)
  q=ifelse(is.vector(lambda), 1,  dim(lambda)[2] )
  psilambda <- 1/psi*lambda 
  M <- diag(q) + crossprod(lambda, psilambda) 
  invM <- chol2inv(chol(M))             
  invsig <- diag(1/psi)- tcrossprod( tcrossprod(psilambda, invM), psilambda) 
  #denoms <- diag(q)- crossprod(lambda, crossprod(invsig, lambda) )
  logdetsig <-  sum(log(psi)) - log(det(invM)) 
  logconst  <- -0.5 * p * log(2*pi) - 0.5 * logdetsig
  X <- sweep(X, 2, mu)
  quadratic_term <-  rowSums(tcrossprod( X , invsig) * X)
  #####
  lpdf<-  logconst  - 0.5 * quadratic_term
  return(lpdf) # a vector
}
################################################################

############################
# gmm.fa2: this uses kmeans initialization. 
############################
library(fad)
gmm.fa.kq <- function(X, K, q , tol = 1e-6, maxiter = 500, nstart= 50){
  X <- as.matrix(X)
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
  #clusters <- kmeans_result$cluster
  #n.clusters <- table(clusters)
  #n.clusters <- as.vector(n.clusters)
  
  lambda <- list() #array(dim = c(p, q, K))
  psi    <- matrix(nrow =K, ncol= p  )  
  for (k in 1:K){ #cov2cor here might not be necessary.
    sigmaR <- cov2cor(sigma[,,k])
    vars <- diag(sigma[,,k])
    start <- 1/ vars
    fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q[k], start= start, maxit = 500, tol = tol)) 
    sd = sqrt(vars)
    lambda[[k]] <-  lmd <-  fa$loadings 
    lambda[[k]]  <- sd * lmd
    psi[k,] <- ps <-  fa$uniquenesses
    psi[k,] <-  sd^2 * ps
  }
  
  llh <- - Inf
  gamma <-  matrix(0, n, K)
  for (t in 1:maxiter){
    # E - Step
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[[k]], psi[k,]))
    logp  <- sweep(logp, 2, log(w), "+")
    maxvals <- apply(logp, 1, max)
    logp  <- sweep(logp, 1, maxvals, "-")
    wprob <- exp(logp)
    sums  <- rowSums(wprob)
    gamma <- sweep(wprob, 1, sums , "/")
    # CM - Step 1: compute mu and weights.
    gsum <- colSums(gamma)
    w <- gsum / n
    matsum <- matrix(0, p, p)
    for (k in 1:K){
      mu[k,] <- colSums(gamma[,k]*X) / gsum[k] 
      # CM - Step 2: compute lambdas and psis
      matsum <- matrix(0, p, p)
      for (i in 1:n){
        covmatt <-  as.matrix(X[i,]- mu[k,])
        prods <-  gamma[i,k] * tcrossprod(covmatt)
        matsum <- matsum + prods
      }
      sigma[,,k] <- matsum /  gsum[k]
    }
    

    clusters <- apply(gamma, 1, function(x) which.max(x) )
    n.clusters <- table(clusters)
    n.clusters <- as.vector(n.clusters)
    
    for (k in 1:K){
      sigmaR <- cov2cor(sigma[,,k])
      vars <- diag(sigma[,,k])
      start <- 1/ vars
      fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q[k], start= start, maxit = 500, tol = tol)) 
      sd = sqrt(vars)
      lambda[[k]] <-  lam <-  fa$loadings 
      lambda[[k]] <- sd * lam # in covariance scale
      psi[k,] <- ps <-  fa$uniquenesses   # We can improve convergence speed by using the updated lambda to compute psi. 
      #First update sigma[,,k]= lambda[,,k]%*%lambda[,,k] + psi[k,] for each k. 
      psi[k,] <- sd^2*ps
    }
    llh0 <- llh
    #####################
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[[k]], psi[k,]))
    logp  <- sweep(logp, 2, log(w), "+")
    maxvals <- apply(logp, 1, max)
    logp  <- sweep(logp, 1, maxvals, "-")
    wprob <- exp(logp)
    sums  <- rowSums(wprob)
    llh   <- sum(maxvals, log(sums))
    # stopping criterion. 
    if  ( abs(llh-llh0) < tol * abs(llh)){
      break
    }
  }
  gamma <- sweep(wprob, 1, sums , "/")
  clusters <- apply(gamma, 1, function(x) which.max(x) )
  iter <- t
  clusters <- apply(gamma, 1, function(x) which.max(x) )
  iter <- t
  
  #nparam <-  (K - 1) + 2 * K * p + K * (p * q - q * (q - 1)/2)
  nparam <-  (K - 1) + K * p + sum (p * q + p - q * (q - 1)/2) 
  BIC <- -2*llh +  nparam * log(n) 
  
  analyser = list(clusters = clusters,  weights = w, 
                  lambda= lambda, psi=psi, sigma = sigma, niter =t,
                  means = mu, loglik= llh, BIC = BIC, #loglik.gmm= gm.mod$loglik, 
                  nfactors = q, converged = t < maxiter, gamma= gamma
  )
  return(analyser)
}  



############################
# gmm.fa3: this uses random initialization.  
############################
library(fad)
gmm.fa.q <- function(X, K, q , tol =1e-6, maxiter = 500, nstart= 100, init=NULL){
  X <- as.matrix(X)
  #q= c(q[1], ..., q[K]) 
  n= nrow(X)
  p= ncol(X)
  if (is.null(init)) {
    set.seed(12321)
  } else {
    set.seed(init)
  }
  
  mu <- matrix(rnorm(K*p),K )
  
  lambda <- list()     # array(dim = c(p, q, K) )
  psi    <- matrix(nrow = K, ncol= p  )  
  sigma  <- array(dim = c(p, p, K) )
  for (k in 1:K){    
    lambda[[k]] <- matrix(rnorm(p*q[k]), p, q[k])
    psi[k,]     <- runif(p, 0.3, 0.9)
    sigma[,,k]  <-  tcrossprod(lambda[[k]]) + diag(psi[k,]) 
  }
  pr<- abs(rnorm(K)) 
  w<- rep(1/K, K)                   # pr/ sum(pr) 
  
  #initials <- list(mean = mu, Sigma =sigma, weights = w , lambda =lambda, psi = psi)
  
  llh <- - Inf
  gamma <-  matrix(0, n, K)
  for (t in 1:maxiter){
    # E - Step
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[[k]], psi[k,]))
    logp  <- sweep(logp, 2, log(w), "+")
    maxvals <- apply(logp, 1, max)
    logp  <- sweep(logp, 1, maxvals, "-")
    wprob <- exp(logp)
    sums  <- rowSums(wprob)
    gamma <- sweep(wprob, 1, sums , "/")
    
    # CM - Step 1
    gsum <- colSums(gamma)
    w <- gsum / n
    matsum <- matrix(0, p, p)
    for (k in 1:K){
      mu[k,] <- colSums(gamma[,k]*X) / gsum[k]
      matsum <- matrix(0, p, p)
      for (i in 1:n){
        covmatt <-  as.matrix(X[i,]- mu[k,])
        prods <-  gamma[i,k] * tcrossprod(covmatt)
        matsum <- matsum + prods
      }
      sigma[,,k] <- matsum /  gsum[k]
    }
    # CM - Step 2: compute lambdas and psis
    for (k in 1:K){
      sigmaR <- cov2cor(sigma[,,k])
      vars <- diag(sigma[,,k])
      start <- 1 / vars
      fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q[k], start= start, maxit = 500, tol = tol)) 
      sd = sqrt(vars)
      lambda[[k]] <-  lam <-  fa$loadings 
      lambda[[k]] <- sd * lam # in covariance scale
      psi[k,] <- ps <-  fa$uniquenesses   # We can improve convergence speed by using the updated lambda to compute psi. 
      #First update sigma[,,k]= lambda[,,k]%*%lambda[,,k] + psi[k,] for each k. 
      psi[k,] <- sd^2*ps
    }
    llh0 <- llh
    #####################
    logp  <- sapply(1:K, function(k) logpdf(X, mu[k,], lambda[[k]], psi[k,]))
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
  
  #nparam <- K-1 + K * p + K * p* (p+1) / 2 + K*(p*q + p) # number of parameters 
  #nparam <- (K - 1) + K * p + p + K * (p * q - q * (q - 1)/2) 
  nparam <-  (K - 1) + K * p + sum (p * q + p - q * (q - 1)/2)
  BIC <- -2*llh +  nparam * log(n) 
  
  analyser = list(clusters = clusters,  weights = w, 
                  lambda= lambda, psi=psi, sigma = sigma, niter =t,
                  means = mu, loglik= llh, BIC = BIC, #loglik.gmm= gm.mod$loglik, 
                  nfactors = q, converged = t < maxiter, gamma= gamma
  )
  return(analyser)
}


##################
# GMM wrapper.new: GMM with multiple initials and inner gmm.fa2.new; (for comparison purposes only)
##################
gmm.fad.q <- function(Y, K, q , tol = 1e-6, maxiter = 500, nstart= 20, init = NULL){
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  #set.seed(123456)
  llhs  <- numeric( (nstart) )  
  res0  <- try({gmm.fa.kq(Y, K, q, tol =tol,  maxiter = 1000 ) }, TRUE)
  llhs0 <- tryCatch({res0$loglik} , error= function(e) -Inf ) 
  
  if (is.null(init)) {
    set.seed(123456)
  } else {
    set.seed(init)
  }  
  seed <- runif(nstart)*10^8
  for (j in 1:nstart){
    res <- try({ gmm.fa.q(Y, K, q, tol =tol, init= seed[j], maxiter = 5) }, TRUE)
    llhs[j] <- tryCatch({res$loglik},  error= function(e) -Inf )     
  }
  llhs1 <- c(llhs, llhs0)
  if (sum(llhs1==-Inf)== (nstart + 1)) stop("All initial values and kmeans encountered errors. 
                      Number of factors or number of clusters might have 
                        been chosen too high compared to sample size.\n
                        Also, try standardizing data. \n")
  irank <- order(llhs, decreasing = TRUE)[1:min(5, max(1,floor(nstart/4)))]
  
  for (j in irank){
  
    res2 <- try({
      gmm.fa.q(Y, K, q,  tol =tol, init= seed[j], maxiter = maxiter) 
    }, TRUE)
    llhs[j] <- tryCatch({res2$loglik},  error= function(e) -Inf)     
  }
  llhs1[1:nstart] <- llhs
  bestindex<- which.max(llhs1)
  if (bestindex== (nstart + 1) ){
    return( res0 )
  } else {
    res.f <-  gmm.fa.q(Y, K, q, tol = tol, init = seed[bestindex],  maxiter = maxiter)
    res.f$init = seed[bestindex]
    return( res.f )
  }
}
#############################


