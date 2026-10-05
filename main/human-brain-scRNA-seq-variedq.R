source("~/Desktop/gmfad paper/gmfad_qq.R")
## ----------------------------------------------------------------------------------
library(doParallel)
library(MixSim) 
library(aricode)
library(radviz3d)
library(cluster)

## load dataset
load("~/Desktop/gmfad paper/darmanis-adult-neuron-astrocyte-200HVG-GDT.rda")

type <- darmanis_adult_2class_df$cell_type
true_cl <- as.integer(type)
table(true_cl)
X <-as.matrix(darmanis_adult_2class_df[, -1])
dim(X)

cc <- detectCores()
cl <- makePSOCKcluster(cc - 1)
registerDoParallel(cl)

set.seed(1233345)

kvals=c(2,3,4)
q1_vals <- c(15,16,17)
q2_vals <- c(12,13,14) 

totals <- length(q1_vals) * length(q2_vals) * length(kvals)

results <- data.frame(
  K = rep(kvals[1],totals),
  q1   = numeric(totals),
  q2   = numeric(totals),
  AIC = rep(NA, totals),
  
  BIC = rep(NA, totals),
  silhouette = rep(NA, totals),
  
  ARI = rep(NA, totals),
  NMI = rep(NA, totals),
  FMI = rep(NA, totals),
  loglik = numeric(totals),
  runtime = numeric(totals),
  convergence= rep(NA, totals),
  niter = numeric(totals)
) 

tol <- 1e-6

idx <- 1

for (K in kvals) {
  for (q1 in q1_vals) {
    for (q2 in q2_vals) {
    
    t1 <- proc.time()
    qvec = c(q1,q2)

    gmfads <- tryCatch({
      gmm.fad.q(X, K, qvec, maxiter = 1000, tol =tol, nstart=20)
    }, error = function(e) NA)
    
    cpu_time <- (proc.time() - t1)[3]
    
    saveRDS(gmfads, file = paste0("scRNA_Tasic_gmmfad_outputParams_output_BIC_K_",K,"_q1_",q1,"_q2_",q2,".rds"))
    
    # Extract outputs
    est_cl <- tryCatch({gmfads$clusters}, error = function(e) NA)
    ari <- tryCatch({ RandIndex(true_cl, est_cl)$AR }, error = function(e) NA)
    bic <- tryCatch({ gmfads$BIC }, error = function(e) NA)
    aic <- tryCatch({ gmfads$AIC }, error = function(e) NA)
    niter_val <- tryCatch({ gmfads$niter }, error = function(e) NA)
    loglik_val <- tryCatch({ gmfads$loglik }, error = function(e) NA)
    # Optional metrics (if functions exist)
    nmi <- tryCatch({ NMI(true_cl, est_cl) }, error = function(e) NA)
    fmi <- tryCatch({ RandIndex(true_cl, est_cl)$F }, error = function(e) NA)
    sil_vals <- tryCatch({ silhouette(est_cl, dist(X)) }, error = function(e) NA)
    sil_mean <- tryCatch({ mean(sil_vals[, 3]) }, error = function(e) NA)
    convergedInd <- tryCatch({ gmfads$converged }, error = function(e) NA)
    # Store results
    results[idx, ] <- c(
      K, q1, q2, aic, bic, sil_mean, 
      ari, nmi, fmi, loglik_val, cpu_time, convergedInd ,niter_val
    )
    
    idx <- idx + 1
  }
}
}