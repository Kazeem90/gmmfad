# The Lymphoma dataset
source("~/Desktop/Dimensionality reduction research/codes/em4gmm.R")
#install.packages("spls")
library(spls)
library(doParallel)
library(MixSim) # for ARI  
data(lymphoma)

cc <- detectCores()
cl <- makePSOCKcluster(cc - 1)
registerDoParallel(cl)

 #str(lymphoma)
#lymphoma$y
lymph <- lymphoma$x
true_cl <- lymphoma$y+1

#q=2 
#gmfads <- gmm.fad2(x, 3, q)
#est_cl <- gmfads$clusters
#ari1 <- RandIndex(true_cl, est_cl); ari1
q1= 18:5
q2= 9:5
q3= 9:5
K = 3
q_grid <- as.matrix(expand.grid(q1=q1, q2=q2, q3=q3))

BICs<-runtimes<- aris<-loglik<- tol <- niter<- numeric(nrow(q_grid) )
K=3 # fixed
tol = 1e-8
tols = rep(1e-8,nrow(q_grid))
print(Sys.time())
for (j in 1:nrow(q_grid) ){
  t1 <- proc.time()
  qvec <- as.vector(q_grid[j,])
  gmfads <- tryCatch({
    gmm.fad.q(lymph, K, qvec, tol =tol, nstart=20)
    }, error= function(e) NA)
  runtimes[j]<- (proc.time()-t1)[3]
  est_cl     <- tryCatch({gmfads$clusters  }, error= function(e) NA)
  aris[j]    <- tryCatch({RandIndex(true_cl, est_cl)$AR  }, error= function(e) NA)
  BICs[j]    <- tryCatch({gmfads$BIC   },     error= function(e) NA)
  niter[j]   <- tryCatch({gmfads$niter  },    error= function(e) NA)
  loglik[j]  <- tryCatch({gmfads$loglik  },   error= function(e) NA)
  #cat("ARI for q=",qvec," is ", aris[j],"\n" )
  #cat("time for q=",qvec," is ", runtimes[j],"\n" )
  #cat("BIC for q=",qvec," is ", BICs[j],"\n" )
  #cat("loglikelihood for q=",qvec," is ",gmfads$loglik ,"\n" ) 
  #cat("no. of iterations for q=",qvec," is ",gmfads$niter,"\n" ) 
  #cat("cluster table is", as.vector(table(gmfads$clusters)))
  
  print(Sys.time())
  
  result_lymphoma <- data.frame(ARI = aris[j], runtimes = runtimes[j], BIC = BICs[j], 
                                niter = niter[j], loglik = loglik[j], tolerance= tols[j], 
                                q1=qvec[1], q2=qvec[2], q3=qvec[3])
  #write.csv(result_lymphomatest, "lymphoma_data_clustering_results_new.csv" )
  print(result_lymphoma)
  
  write.table( result_lymphoma, file = "lymphoma_data_clustering_results_new2.csv" , append = TRUE,  sep = ",", 
               col.names = !file.exists("lymphoma_data_clustering_results_new2.csv"), row.names = FALSE )
}

result_lymphoma <- data.frame(ARI = aris, runtimes = runtimes, BIC = BICs, 
                              niter = niter, loglik = loglik, tolerance= tols)
print(result_lymphoma)
#result_lymphomatest <- data.frame(ARI = aris, runtimes = runtimes, BIC = BICs)

# dist_matrix <- as.dist(data_matrix)

# Perform MDS
# mds_result <- cmdscale(dist_matrix, k = 2) # k is the number of dimensions

# Display the result
#print("MDS Result:")
#print(mds_result)
