###################################################################################
###########################################################
##### BIC Simulations with gmm factor analyzers: gmm.fad vs EMMIXmfa: Optimal K detection rate. 
###########################################################
###################################################################################
## simulation with gmm factor analysers
source("~/Desktop/Dimensionality reduction research/codes/em4gmm.R")
source("~/Desktop/Dimensionality reduction research/codes/gmfad_qq.R")

library(doParallel)
library(EMMIXmfa)
library(MASS)
library(mclust)
library(MixSim)
library(clue)
library(ggplot2)
library(tidyr)
library(dplyr)
# function for determining level of separation of dataset.
get_omega <- # function to obtain alpha value given a specified omega
  function(w_true,mu_true, sigma_true,omega = 0.001, p,
           tolerance= min(omega*1e-2, 1e-4*p), 
           alpha= 1, alpha_steps =1e-4 ) {
    if (p >20) alpha_steps = 5*p*alpha_steps 
    omega_hat <- MixSim::overlapGOM(Pi = w_true, Mu = mu_true, S = alpha * sigma_true)
    i=0
    while ( abs(omega_hat - omega) > tolerance && alpha>0 && i<1e5){
      alpha = ifelse( (omega_hat - omega) > 0, (alpha - alpha_steps), (alpha + alpha_steps) ) 
      omega_hat <- MixSim::overlapGOM(Pi = w_true, Mu = mu_true, S = alpha * sigma_true); 
      i=i+1
      #if (i %% 100==0)    print(c(alpha= alpha, omega= omega_hat, i= i)) 
      if ( omega_hat < omega | i <=1) {
        alpha = alpha + 1
        omega_hat = 1
      }
    }
    if (alpha <= 0) stop("alpha cannot be negative. Try smaller alpha_steps")
    return( list(omega_hat=omega_hat, alpha = alpha, niter=i))
  }

cc <- detectCores()
cl <- makePSOCKcluster(cc - 1)
registerDoParallel(cl)
n_sim =  100
omegavals = c(0.001, 0.005, 0.01) 
ptm <-proc.time(); print(Sys.time())
seedling = 13579
set.seed(seedling)
seeds <- runif(n_sim) * 10^8
inits <- runif(n_sim) * 10^8 # for tracking errors
n <- 300 
kvals =  c(2) 

pp = c(10) #, 50, 100, 200) 

for (K in kvals){
  print(Sys.time())
  G <- 2*K
  q<- c(2,3)
  correctrates1 <- matrix(0,length(omegavals), length(pp))
  failrates1 <- failrates2 <-  matrix(0,length(omegavals), length(pp) )
  rownames(correctrates1) <-  rownames(failrates1) <- c("alpha0.001", "alpha0.005", "alpha0.01")
  for (omega in omegavals){
    print(Sys.time())
    
    model_correct <-  modelfail1 <- matrix(0, n_sim, length(pp))
    for (p in pp){
      Z <- array(dim=c(n,p, n_sim))
      Zid <- array(dim= c(n, n_sim))
      # w_true <- matrix(0, n_sim, K) 
      for (i in 1:n_sim){
        if(i %% 20 == 1) cat("gmmfad omega=",omega,", n =", n, ", p =", p,", q =", q, ", K =", K,", it =", i, "\n")
        set.seed(seeds[i])
        prob <- runif(K, 0.2, 0.9)
        w_true = prob/sum(prob)
        psi_true <- matrix( runif(K*p, 0.3, 0.9),   nrow = K, ncol = p)
        mu_true <- matrix( rnorm(K*p), nrow = K, ncol = p)
        lambda_true <- list()  
        sigma_true  <- array(dim = c(p, p, K) )
        for (k in 1:K){    
          lambda_true[[k]] <- matrix(rnorm(p*q[k]), p, q[k])
          sigma_true[,,k]  <- tcrossprod(lambda_true[[k]]) + diag(psi_true[k,]) 
        }
        
        alpha <-get_omega(w_true,mu_true, sigma_true, omega = omega, p=p)$alpha
        Ysim <- simdataset(n = n, Pi = w_true, Mu = mu_true, S= alpha * sigma_true)
        Z[,,i] <- Ysim$X
        Zid[,i] <- Ysim$id
      } # data generation ends here.
      
      cat("gmmfad omega=",omega," data generated for  n =", n, ", p =", p,", q =", q, ", K =", K, "\n")
      for (i in 1:n_sim){
        Y <- Z[,,i]
        #if(i %% 5 == 1)
        cat("gmmfad omega=",omega,", n =", n, ", p =", p,", q =", q, ", K =", K,", it =", i, "\n")
        bic1 <- array(NA, dim = c( (2*q[1]), (2*q[2]), (2*K) ) )
        for (g in 1:G){
          for (qq2 in 1:(2*q[2])){
            for (qq1 in 1:(2*q[1])){
              if (g==1){
                qq = c(qq1) 
              } else if(g==2){
                qq = c(qq1, qq2)
              } else if(g==3){
                qq = c(qq1, qq2, qq1)
              }else{
                qq = c(qq1, qq2, qq1, qq2)
              }
              gmmfa <- tryCatch({  gmm.fad.q(Y, K=g, q=qq, maxiter = 500, nstart= 20, init= seeds[i] )
              }, error = function(e) NA)
              bic1[qq1, qq2, g]  <-  tryCatch({ gmmfa$BIC }, error = function(e) NA)
            }
          }
        }
        best_params <-  which(bic1 == min(bic1, na.rm = T), arr.ind = TRUE)
        modelfail1[i, which(pp==p)] <- sum(is.na(bic1)) != 0 
        model_correct[i, which(pp==p)] <- any(apply(best_params, 1, function(row) { setequal(row[1:K],q)  && row[(K+1)]== K }))
        cat("correctness for current iteration  is", model_correct[i, which(pp==p)], "\n")
      }
      correctrates1[which(omegavals==omega), which(pp==p) ] <- mean(model_correct[ , which(pp==p)], na.rm = T)  
      failrates1[which(omegavals==omega), which(pp==p)] <- mean(modelfail1[ , which(pp==p)], na.rm = T )
      cat("correctrate for  gmfad omega =",omega,"p=",p ,", q=",q ," is ",correctrates1[which(omegavals==omega), which(pp==p)],"\n")
      cat("failrate for gmfad omega =",omega,"  p=",p ,", q=",q ,"K=",K," is ",failrates1[which(omegavals==omega), which(pp==p)], "\n")
      #print( best_params )
      print(Sys.time())
    }
    #write.csv(data.frame(best_params), paste0("gmfad_bestparams_omega_",omega, "_K_",K, "_q_",q, "_n_",n,".csv"))
  }
  cat("correctrates for gmfad   p=",p ,", q=",q ,"K=",K, " is ", t(correctrates1))
  cat("failrates for    gmfad   p=",p ,", q=",q ,"K=",K, " is ", t(failrates1))
  print(Sys.time())
  write.csv(data.frame(t(correctrates1)), paste0("gmfad2_q24_BICcorrectrates_omega_K_",K,"_n_",n, ".csv"))
  write.csv(data.frame(t(failrates1)),    paste0("gmfad2_q24_BICfailrates_omega_K_",K,"_n_",n, ".csv"))
}