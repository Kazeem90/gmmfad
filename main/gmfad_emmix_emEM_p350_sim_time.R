source("~/Desktop/Dimensionality reduction research/codes/em4gmm.R")

library(clue)
library(MixSim)
############# funnctions
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
        omega_hat =1
      }
    }
    if (alpha <= 0) stop("alpha cannot be negative. Try smaller alpha_steps")
    return( list(omega_hat=omega_hat, alpha = alpha, niter=i))
  }
recoloring_gmfad  <- function(estimated_clusters, true_clusters, mu, lambda, psi){
  contingency_table <- table(true_clusters, estimated_clusters)
  assignment <- solve_LSAP(contingency_table, maximum = TRUE)
  mapping <- as.vector(assignment)
  reordered_mu<- mu
  reordered_lambda<-lambda
  reordered_psi <- psi
  for (k in 1:K){
    reordered_mu[k,]<- mu[mapping[k], ]
    reordered_lambda[,,k] <- lambda[, , mapping[k]]
    reordered_psi[k,]<- psi[mapping[k], ]
  }
  return(list(mu= reordered_mu, lambda=reordered_lambda, 
              psi=reordered_psi, mapping=mapping) )
}  
recoloring_mfa  <- function(estimated_clusters, true_clusters, mu, lambda, psi){
  contingency_table <- table(true_clusters, estimated_clusters)
  assignment <- solve_LSAP(contingency_table, maximum = TRUE)
  mapping <- as.vector(assignment)
  reordered_mu<- mu
  reordered_lambda<-lambda
  reordered_psi <- psi
  for (k in 1:K){
    reordered_mu[,k]      <- mu[, mapping[k]]
    reordered_lambda[,,k] <- lambda[,,mapping[k]]
    reordered_psi[,,k]    <- psi[, , mapping[k]]
  }
  return(list(mu= reordered_mu, lambda=reordered_lambda, 
              psi=reordered_psi, mapping=mapping) )
}

########################################################################################################
##############################################################
### Simulations with gmm factor analysers: gmm.fad vs EMMIXmfa
library(doParallel)
library(EMMIXmfa)
library(MASS)
library(mclust)
library(MixSim)
library(clue)
library(ggplot2)
library(tidyr)
library(dplyr)
cc <- detectCores()
cl <- makePSOCKcluster(cc - 1)
registerDoParallel(cl)
#omegavals <- c(0.001, 0.005, 0.01)
#for (omega in omegavals){
  n_sim = 100
  #omega = 0.001 
  omega = NA
  ptm <-proc.time(); Sys.time()
  seedling= 13579
  set.seed(seedling)
  seeds <- runif(n_sim) * 10^8
  inits <- runif(n_sim) * 10^8 # for tracking errors
  n = 150 
  # K = 2 # 
  kvals = c(2,3)#, 3)
  qvals = c(3,2) #2, 3)
  pp = c(150) #, 50, 100, 200)     #, 1000)
  timelist.k  <- list()
  ARIlist.k   <- list()
  frobs.mu.k  <- list()
  frobs.lambda.k<- list()
  frobs.psi.k <- list()
  
  
  for (K in kvals){
    
    timelist.q <- list()
    ARIlist.q   <- list()
    frobs.mu.q <- list()
    frobs.lambda.q <- list()
    frobs.psi.q <- list()
    
    
    for (q in qvals){
      
      timelist.p <- list()
      ARIlist.p   <- list()
      frobs.mu.p<- list()
      frobs.lambda.p<- list()
      frobs.psi.p <- list()
    
      
      for (p in pp){
        seeds <- runif(n_sim) * 10^8
        inits <- runif(n_sim) * 10^8
        Z   <- array(dim=c(n,p, n_sim))
        Zid <- array(dim= c(n, n_sim))
        t.gmfad <- t.emmix  <- rel.speedup <-  numeric(n_sim)
        ar.gmfad <- ar.emmix <- numeric(n_sim)
        frob.mu.gmfad <- frob.mu.emmix<-  matrix(0, n_sim, K)
        frob.lambda.gmfad <- frob.lambda.emmix<-  matrix(0, n_sim, K)
        frob.psi.gmfad <- frob.psi.emmix <- matrix(0, n_sim, K)
        for (i in 1:n_sim){
          if(i %% 10 == 1) cat("Data simulation nsim =",n_sim,", n =", n, 
                               ", p =", p,", q =", q, ", K =", K,", omega=",omega,", it =", i, "\n")
          
          set.seed(seeds[i])
          prob <-runif(K, 0.2, 0.9)
          w_true = prob / sum(prob)
          psi_true <- matrix(runif(K*p, 0.3, 0.9),  nrow = K, ncol = p)  
          mu_true <- matrix( rnorm(K*p), nrow = K, ncol = p)
          lambda_true <- array(rnorm(p*q*K), dim = c(p,q,K) )
          sigma_true <- array( dim = c(p,p,K) )
          for (k in 1:K){
            sigma_true[,,k] <- tcrossprod(lambda_true[,,k], lambda_true[,,k]) + diag(psi_true[k, ])
          }
          # simulate mixed data using the above parameters
          alpha <-1 #get_omega(w_true,mu_true, sigma_true, omega = omega, p=p)$alpha
          lambda_true = sqrt(alpha) * lambda_true
          psi_true = alpha * psi_true
          Xsim <- simdataset(n = n, Pi = w_true, Mu = mu_true, S= alpha * sigma_true)
          Y <- Xsim$X 
          cluster_true = Xsim$id
          print(Sys.time())
          
          #if (i %% 5 == 1) 
          cat("gmmfad vs EMMIXmfa nsim=",n_sim,", n =", n, ", p =", p,", q =", q, ", K =", K,",omega=",omega, ",it =", i, "\n")
          prc <- proc.time()
          gmfad <- gmm.fad(Y, K=K, q=q, maxiter = 500, nstart= 20, tol = 1e-6, init= seeds[i] ,  innerNStart= 3)
          t.gmfad[i] <- (proc.time() - prc)[3]
          cluster_est <- gmfad$clusters
          ar.gmfad[i] <- MixSim::RandIndex(cluster_true, cluster_est)$AR
          # reordering estimated parameters before computing average Frobenius distance.  
          reordered_params <- recoloring_gmfad(cluster_est, cluster_true, gmfad$means, gmfad$lambda,gmfad$psi)
          mu_est<- reordered_params$mu
          lambda_est<- reordered_params$lambda
          psi_est <- reordered_params$psi
          
          for (k in 1:K){
            frob.mu.gmfad[i,k ] <-  
              norm(as.matrix(mu_est[k,]-mu_true[k,]), type ="F" ) / norm(as.matrix(mu_true[k,]), type ="F" )
            frob.lambda.gmfad[i,k ] <-  
              norm((tcrossprod(lambda_est[,,k]) - tcrossprod(lambda_true[,,k])), type ="F" ) / norm(tcrossprod(lambda_true[,,k]), type ="F")
            frob.psi.gmfad[i,k ] <-  
              norm(as.matrix(psi_est[k,]-psi_true[k,]), type ="F" ) / norm(as.matrix(psi_true[k,]), type ="F")
          }
          ### emmixmfa
          prc <- proc.time()
          emmix <- emmix_emEM(Y, K=K, q=q, tol = 1e-6, maxiter = 500, nstart= 20, init = seeds[i],  innerNStart= 3)#D_type = 'unique', sigma_type='unique', tol = 1e-6, 
          t.emmix[i] <- (proc.time() - prc)[3]
          rel.speedup[i] <- t.emmix[i] / t.gmfad[i]
          
          if (!is.null(emmix) ){  cluster_est <- emmix$clust
          } else{  cluster_est <- NA 
          }
          ar.emmix[i] <- tryCatch({ MixSim::RandIndex(cluster_true, cluster_est)$AR  },  error= function(e) NA) 
          # reordering estimated parameters before computing average Frobenius distance.  
          reordered_params <- tryCatch({ recoloring_mfa(cluster_est, cluster_true, emmix$mu, emmix$B,emmix$D) 
          }, error= function(e) NA) ###
          mu_est<- tryCatch({  reordered_params$mu },  error= function(e) NA) 
          lambda_est<- tryCatch({ reordered_params$lambda },  error= function(e) NA) 
          psi_est <- tryCatch({ reordered_params$psi },  error= function(e) NA) 
          
          for (k in 1:K){
            frob.mu.emmix[i,k ] <-  tryCatch({ 
              norm(as.matrix(t(mu_est)[k,]-mu_true[k,]), type ="F" ) / norm(as.matrix(mu_true[k,]), type ="F" )
            },  error= function(e) NA) 
            frob.lambda.emmix[i,k ] <-  tryCatch({ 
              norm((tcrossprod(lambda_est[,,k]) -tcrossprod(lambda_true[,,k])), type ="F" ) / norm(tcrossprod(lambda_est[,,k]) , type ="F" )
            },  error= function(e) NA) 
            frob.psi.emmix[i,k ] <-  tryCatch({ 
              norm(as.matrix(diag(psi_est[,,k])-psi_true[k,]), type ="F" ) / norm(as.matrix(psi_true[k,]), type ="F" )
            },  error= function(e) NA) 
          }
        } 
        
        times       <- data.frame(gmfad=t.gmfad, emmix = t.emmix, speedup.rate= rel.speedup )
        ARI         <- data.frame(gmfad=ar.gmfad, emmix = ar.emmix)
        frob.mu     <- as.data.frame(cbind(frob.mu.gmfad,    frob.mu.emmix))
        frob.lambda <- as.data.frame(cbind(frob.lambda.gmfad,frob.lambda.emmix))
        frob.psi    <- as.data.frame(cbind(frob.psi.gmfad,   frob.psi.emmix))
        #
        timelist.p[[which(pp==p)]] <- times # update to timelist.p
        ARIlist.p[[ which(pp==p)]] <- ARI
        frobs.mu.p[[which(pp==p)]]<- frob.mu
        frobs.lambda.p[[which(pp==p)]]<- frob.lambda
        frobs.psi.p[[which(pp==p)]]<- frob.psi
        
        write.csv(times, paste0("gmmfa_emmix_time_nsim_",n_sim,"_n_",n,"_p_",p,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
        write.csv(ARI, paste0("gmmfa_emmix_ARI_nsim_",n_sim,"_n_",n,"_p_",p,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
        write.csv(as.data.frame(seeds), paste0("gmmfa_emmix_seeds_nsim_",n_sim,"_n_",n,"_p_",p,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
        write.csv(frob.mu, paste0("gmmfa_emmix_Frobenius_mu_nsim_",n_sim,"_n_",n,"_p_",p,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
        write.csv(frob.lambda, paste0("gmmfa_emmix_Frobenius_lambda_nsim_",n_sim,"_n_",n,"_p_",p,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
        write.csv(frob.psi, paste0("gmmfa_emmix_Frobenius_psi_nsim_",n_sim,"_n_",n,"_p_",p,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
      } # loop out of p
      print(Sys.time())
      
      for (j in 1:length(pp)){
        timelist.p[[j]]$p      = rep(pp[j], nrow(timelist.p[[1]]))
        ARIlist.p[[j]]$p       = rep(pp[j], nrow(timelist.p[[1]]))
        frobs.mu.p[[j]]$p      = rep(pp[j], nrow(timelist.p[[1]]))# replace n_sim with length(ARIlist1[[j]]$p)
        frobs.lambda.p[[j]]$p  = rep(pp[j], nrow(timelist.p[[1]]))
        frobs.psi.p[[j]]$p     = rep(pp[j], nrow(timelist.p[[1]]))
      }
      # merge
      timedf.p         <-  do.call(rbind,  timelist.p) 
      ARIdf.p          <-  do.call(rbind,  ARIlist.p) 
      frobs.mu.df.p    <-  do.call(rbind,  frobs.mu.p)
      frobs.lambda.df.p<-  do.call(rbind,  frobs.lambda.p)
      frobs.psi.df.p   <-  do.call(rbind,  frobs.psi.p)
      #
      write.csv(timedf.p, paste0("gmmfa_emmix_all_time_nsim_",n_sim,"_n_",n,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
      write.csv(ARIdf.p, paste0("gmmfa_emmix_all_ARI_nsim_",n_sim,"_n_",n,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
      write.csv(frobs.mu.df.p, paste0("gmmfa_emmix_all_frobs_mu_nsim_",n_sim,"_n_",n,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
      write.csv(frobs.lambda.df.p, paste0("gmmfa_emmix_all_frobs_lambda_nsim_",n_sim,"_n_",n,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
      write.csv(frobs.psi.df.p, paste0("gmmfa_emmix_all_frobs_psi_nsim_",n_sim,"_n_",n,"_K_",K,"_q_",q,"_omega_",omega,"_seed_",seedling,".csv"))
      
      timelist.q[[which(qvals==q)]] <- timedf.p # update to timelist.p
      ARIlist.q[[which(qvals==q)]]  <- ARIdf.p 
      frobs.mu.q[[which(qvals==q)]] <- frobs.mu.df.p
      frobs.lambda.q[[which(qvals==q)]]<- frobs.lambda.df.p
      frobs.psi.q[[which(qvals==q)]]<- frobs.psi.df.p
      
    }                    # end of q loop. 
    print(Sys.time())
    
    for (j in 1:length(qvals)){
      timelist.q[[j]]$q    =  rep(qvals[j], dim(timelist.q[[j]])[1])
      ARIlist.q[[j]]$q     =  rep(qvals[j], dim(timelist.q[[j]])[1])
      frobs.mu.q[[j]]$q    =  rep(qvals[j], dim(timelist.q[[j]])[1])# replace n_sim with length(ARIlist1[[j]]$p)
      frobs.lambda.q[[j]]$q=  rep(qvals[j],dim(timelist.q[[j]])[1])
      frobs.psi.q[[j]]$q  =  rep(qvals[j],  dim(timelist.q[[j]])[1])
    } 
    #merge
    timedf.q          <- do.call(rbind,  timelist.q) 
    ARIdf.q           <- do.call(rbind,  ARIlist.q) 
    frobs.mu.df.q     <- do.call(rbind,  frobs.mu.q)
    frobs.lambda.df.q <- do.call(rbind,  frobs.lambda.q)
    frobs.psi.df.q    <- do.call(rbind,  frobs.psi.q)
    ###
    timelist.k[[which(kvals==K)]] <- timedf.q 
    ARIlist.k[[which(kvals==K)]]  <- ARIdf.q 
    #
    frobs.mu.df.q$K <- rep(K, nrow(frobs.mu.df.q) ) 
    frobs.lambda.df.q$K <- rep(K, nrow(frobs.mu.df.q) ) 
    frobs.psi.df.q$K <- rep(K, nrow(frobs.mu.df.q) ) 
    ###
    print(Sys.time())
    write.csv(frobs.mu.df.q, paste0("gmmfa_emmix_emEM_allpq_frobs_mu_nsim_",n_sim,"_n_",n,"K_",K,"_all_pq_omega_",omega,"_seed_",seedling,".csv"))
    write.csv(frobs.lambda.df.q, paste0("gmmfa_emmix_emEM_allpq_frobs_lambda_nsim_",n_sim,"_n_",n,"K_",K,"_all_pq_omega_",omega,"_seed_",seedling,".csv"))
    write.csv(frobs.psi.df.q, paste0("gmmfa_emmix_emEM_allpq_frobs_psi_nsim_",n_sim,"_n_",n,"K_",K,"_all_pq_omega_",omega,"_seed_",seedling,".csv"))
  }                    # end of K loop
  ################################
  ################################
  ################################
  for (j in 1:length(kvals)){
    timelist.k[[j]]$K =  rep(kvals[j], dim(timelist.k[[j]])[1])
    ARIlist.k[[j]]$K  =  rep(kvals[j], dim(timelist.k[[j]])[1])
   
  } 
  # merge
  timedf.k         <- do.call(rbind, timelist.k) 
  ARIdf.k          <- do.call(rbind, ARIlist.k) 
  
 
  write.csv(timedf.k, paste0("gmmfa_emmix_emEM_all_time_nsim_",n_sim,"_n_",n,"_all_pqK_omega_",omega,"_seed_", seedling,".csv"))
  write.csv(ARIdf.k, paste0("gmmfa_emmix_all_emEM_ARI_nsim_",n_sim,"_n_",n,"_all_pqK_omega_",omega,"_seed_",seedling,".csv"))

  Sys.time()
  duration <- proc.time() - ptm; duration
  #}
#}
########################################
########################################

