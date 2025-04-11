###################################################################################
###########################################################
### BIC Simulations with gmm factor analyzers: gmm.fad vs EMMIXmfa: Optimal K detection rate. 
###########################################################
## simulation with gmm factor analysers
source("~/Desktop/Dimensionality reduction research/codes/em4gmm.R")

library(doParallel)
library(EMMIXmfa)
library(MASS)
library(mclust)
library(MixSim)
library(clue)
library(ggplot2)
library(tidyr)
library(dplyr)

get_omega_old <- # function to obtain alpha value given a specified omega
  function(w_true,mu_true, sigma_true,omega = 0.001, 
           tolerance= 1e-6, alpha= 2, alpha_steps =1e-4){
    omega_hat <- 1
    while ( (omega_hat- omega) >tolerance && alpha>0){
      alpha <- alpha - alpha_steps 
      omega_hat <- MixSim::overlapGOM(Pi = w_true, Mu = mu_true, S = alpha * sigma_true) #)
    }
    if (alpha<= 0) stop("alpha cannot be negative. Try smaller alpha_steps")
    return( list(omega_hat=omega_hat, alpha = alpha))
  }

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
      #if (i %% 100==0) print(c(alpha= alpha, omega= omega_hat, i= i)) 
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
n_sim = 100
omega = 0.005 # 0.001, 0.01
ptm <-proc.time(); print(Sys.time())
seedling = 13579
set.seed(seedling)
seeds <- runif(n_sim) * 10^8
inits <- runif(n_sim) * 10^8 # for tracking errors
n <- 300 
kvals = c(2, 3)
qvals = c(2, 3)
pp = c(10)
for (K in kvals){
  print(Sys.time())
  G <- 2*K
  correctrates1 <- correctrates2 <- matrix(0,length(qvals), length(pp))
  failrates1<- failrates2 <- matrix(0,length(qvals), length(pp))
  #rownames(correctrates1) <- c("q_2","q_3")
  for (q in qvals){
    print(Sys.time())
    Q <- 2*q
    bestp1 <- bestp2 <- modelfail1 <- modelfail2 <- matrix(0, n_sim, length(pp))# numeric(n_sim) 
    # colnames(bestp1) <- c("p_10","p_50","p_100","p_200")
    # colnames(bestp2) <- c("p_10","p_50","p_100","p_200")
    for (p in pp){
      Z <- array(dim=c(n,p, n_sim))
      Zid <- array(dim= c(n, n_sim))
      ###
      w_true <- matrix(0, n_sim, K) 
      for (i in 1:n_sim){
        if(i %% 20 == 1) cat("gmmfad omega=",omega,", n =", n, ", p =", p,", q =", q, ", K =", K,", it =", i, "\n")
        set.seed(seeds[i])
        #prob = abs(rnorm(K))
        prob <-runif(K, 0.2, 0.9)
        w_true = prob/sum(prob)
        psi_true <- matrix( runif(K*p, 0.3, 0.9),  nrow = K, ncol = p) 
        mu_true <- matrix( rnorm(K*p), nrow = K, ncol = p)
        lambda_true <- array(rnorm(p*q*K), dim = c(p,q,K) )
        sigma_true <- array( dim = c(p,p,K) )
        for (k in 1:K){
          sigma_true[,,k] <- tcrossprod(lambda_true[,,k], lambda_true[,,k]) + diag(psi_true[k, ])
        }
        alpha <-get_omega(w_true,mu_true, sigma_true, omega = omega, p=p)$alpha
        Ysim <- simdataset(n = n, Pi = w_true, Mu = mu_true, S= alpha * sigma_true)
        Z[,,i] <- Ysim$X
        Zid[,i] <- Ysim$id
      } #data generation ends here.
      ###
      cat("gmmfad omega=",omega," data generated for n =", n, ", p =", p,", q =", q, ", K =", K, "\n")
      for (i in 1:n_sim){
        Y <- Z[,,i]
        if(i %% 5 == 1) cat("gmmfad omega=",omega,", n =", n, ", p =", p,", q =", q, ", K =", K,", it =", i, "\n")
        bic1 <- matrix(nrow= G, ncol = Q)
        for (g in 1:G){
          for (q1 in 1:Q){
            gmmfa <- tryCatch({ gmm.fad(Y, K=g, q=q1, maxiter = 500, nstart= 20, init= seeds[i] )
            }, error = function(e) NA)
            bic1[g,q1] <- tryCatch({gmmfa$BIC }, error = function(e) NA)
          }
        }
        bestp1[i, which(pp==p)] <- which.min(bic1)
        modelfail1[i, which(pp==p)] <- sum(is.na(bic1))!=0 
      }
      correctrates1[which(qvals==q), which(pp==p)]<- mean(bestp1[ , which(pp==p)] == K * (2*q-1) , na.rm = T) 
      failrates1[which(qvals==q), which(pp==p)] <- mean(modelfail1[ , which(pp==p)], na.rm = T )
      
      cat("correctrate for gmfad omega =",omega,"p=",p ,", q=",q ," is ",correctrates1[which(qvals==q), which(pp==p)],"\n")
      cat("failrate for gmfad omega =",omega," p=",p ,", q=",q ,"K=",K," is ",failrates1[which(qvals==q), which(pp==p)], "\n")
      
      print(bestp1)
      for (i in 1:n_sim){
        Y <- Z[,,i]
        cat("correctrate for gmfad omega =",omega,"p=",p ,", q=",q ," is ",correctrates1[which(qvals==q), which(pp==p)],"\n")
        cat("emmix nsim=",n_sim,", n =", n, ", p =", p,", q =", q, ", K =", K,", it =", i, "\n")
        bic2 <- matrix(nrow= G, ncol = Q)
        for (g in 1:G){
          for (q1 in 1:Q){
            emmix <- tryCatch({ mfa(Y, g=g, q=q1, D_type = 'unique', sigma_type='unique', tol = 1e-6, conv_measure = 'ratio')
            }, error = function(e) NA)
            bic2[g,q1] <- tryCatch({ ifelse(is.null(emmix$BIC), NA, emmix$BIC)
            }, error = function(e) NA)
          }
        }
        bestp2[i, which(pp==p)] <- which.min(bic2)
        modelfail2[i, which(pp==p)] <- sum(is.na(bic2)) != 0
      }
      print(Sys.time())
      correctrates2[which(qvals==q), which(pp==p)]<- mean(bestp2[ , which(pp==p)] == K * (2*q-1) , na.rm = T) 
      failrates2[which(qvals==q), which(pp==p)] <- mean(modelfail2[ , which(pp==p)], na.rm = T )
      cat("correctrate for gmfad omega =",omega," p=",p ,", q=",q ,"K=",K, " is ",correctrates1[which(qvals==q), which(pp==p)], "\n")
      cat("correctrate for emmmix omega =",omega," p=",p ,", q=",q ,"K=",K," is ",correctrates2[which(qvals==q), which(pp==p)], "\n")
      cat("failrate for gmfad omega =",omega," p=",p ,", q=",q ,"K=",K," is ",failrates1[which(qvals==q), which(pp==p)], "\n")
      cat("failrate for emmix omega =",omega," p=",p ,", q=",q ,"K=",K," is ",failrates2[which(qvals==q), which(pp==p)], "\n")
    }
    write.csv(data.frame(bestp1), paste0("gmfad_bestp1_omega_",omega, "_K_",K, "_q_",q, "_n_",n,".csv"))
    write.csv(data.frame(bestp2), paste0("emmix_bestp2_omega_",omega, "_K_",K, "_q_",q, "_n_",n,".csv"))
  }
  cat("correctrates for gmfad omega =",omega," p=",p ,", q=",q ,"K=",K, " is ", correctrates1)
  cat("correctrates for emmix omega =",omega," p=",p ,", q=",q ,"K=",K, " is ", correctrates2)
  cat("failrates for gmfad omega =",omega," p=",p ,", q=",q ,"K=",K, " is ", failrates1)
  cat("failrates for emmix omega =",omega," p=",p ,", q=",q ,"K=",K, " is ", failrates2)
  print(Sys.time())
  write.csv(data.frame(correctrates1), paste0("gmfad_BICcorrectrates_omega_", omega, "_K_",K,"_n_",n, ".csv"))
  write.csv(data.frame(correctrates2), paste0("emmix_BICcorrectrates_omega_", omega, "_K_",K,"_n_",n, ".csv"))
  write.csv(data.frame(failrates1), paste0("emmix_BICfailrates_omega_", omega, "_K_",K,"_n_",n, ".csv"))
  write.csv(data.frame(failrates2), paste0("emmix_BICfailrates_omega_", omega, "_K_",K,"_n_",n, ".csv"))
}

