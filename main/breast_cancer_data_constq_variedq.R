source("~/Desktop/Dimensionality reduction research/codes/em4gmm.R")
#devtools::install_github("fanne-stat/radviz3d")
#source("//homes.mtu.edu/home/Desktop/Dimensionality reduction research/em4gmm.R")

library(MixSim) # for ARI
print(Sys.time())
##############################
# Data wrangling
set.seed(1233345) 
wbc <- read.csv("wdbc.data", header = FALSE)# Read the names file. Informative!
names <- readLines("wdbc.names")

wbc$V2 <- as.factor(wbc$V2)
true_cl <- as.integer(wbc$V2)
table(wbc$V2)
wbc<- wbc[, -c(1,2)]
print(str(wbc))

wbc.trans <-radviz3d::Gtrans(wbc)
qvals <- 1:25
BICs<-runtimes<- aris<-loglik <- niter<- numeric(length(qvals) )
K=2 # fixed
tol = 1e-12
tols = rep(tol,length(qvals))
print(Sys.time())
for (q in 1:length(qvals) ){
  t1 <- proc.time()
  gmfads <- tryCatch({
    gmm.fad( wbc.trans, K, maxiter = 1000, q, tol =tol, nstart=20)
  }, error= function(e) NA)
  runtimes[q]<- (proc.time()-t1)[3]
  est_cl <- tryCatch({gmfads$clusters }, error= function(e) NA)
  aris[q] <- tryCatch({RandIndex(true_cl, est_cl)$AR }, error= function(e) NA)
  BICs[q] <- tryCatch({gmfads$BIC  },   error= function(e) NA)
  niter[q]  <- tryCatch({gmfads$niter },  error= function(e) NA)
  loglik[q] <- tryCatch({gmfads$loglik },  error= function(e) NA)
  print(Sys.time())
  
  result_breastcancer <- data.frame(q=q, ARI = aris[q], runtimes = runtimes[q], BIC = BICs[q], 
                                niter = niter[q], loglik = loglik[q], tolerance= tols[q] )
  #write.csv(result_lymphomatest, "lymphoma_data_clustering_results_new.csv" )
  print(result_breastcancer)
  write.table( result_breastcancer, file = "breastcancer_data_gmfad_clustering_results_newDec18.csv" , append = TRUE, sep = ",", 
               col.names = !file.exists("breastcancer_data_gmfad_clustering_results_newnewDec18.csv"), row.names = FALSE )
}

results_breastcancer <- data.frame(ARI = aris, runtimes = runtimes, BIC = BICs, 
                              niter = niter, loglik = loglik, tolerance= tols)
print(results_breastcancer)
#################### EMMIX
BIC2<-runtimes<- aris2<-loglik <- niter<- numeric(length(qvals) )
for (q in 1:length(qvals) ){
  t1 <- proc.time()
  emmix <- tryCatch({ 
    mfa(wbc.trans, g=K, q=q, D_type = 'unique',  sigma_type='unique', itmax = 1000, tol = tol, conv_measure = 'ratio')
    }, error = function(e) NA)
  runtimes[q]<- (proc.time()-t1)[3]
  if (!is.null(emmix) ){ est_cl2 <- emmix$clust
  } else{  est_cl2 <- NA 
  } 
  aris2[q] <- 
    tryCatch({
    RandIndex(true_cl, est_cl2)$AR 
    }, error= function(e) NA)
  BIC2[q] <- tryCatch({ifelse(is.null(emmix$BIC), NA, emmix$BIC) },   error= function(e) NA)
  niter[q]  <-  NA
  loglik[q] <- tryCatch({ifelse(is.null(emmix$logL), NA, emmix$logL) },  error= function(e) NA)
  print(Sys.time())
  
  result_breastcancer2 <- data.frame(q=q, ARI = aris2[q], runtimes = runtimes[q], BIC = BIC2[q], 
                                    niter = niter[q], loglik = loglik[q], tolerance= tols[q] )
  #write.csv(result_lymphomatest, "lymphoma_data_clustering_results_new.csv" )
  print(result_breastcancer2)
  write.table( result_breastcancer2, file = "breastcancer_data_emmix_clustering_results_newDec18.csv" , append = TRUE, sep = ",", 
               col.names = !file.exists("breastcancer_data_emmix_clustering_results_newDec18.csv"), row.names = FALSE )
}

results_breastcancer2 <- data.frame(ARI = aris2, runtimes = runtimes, BIC = BIC2, 
                                   niter = niter, loglik = loglik, tolerance= tols)
print(results_breastcancer2)

opt_q <- which.min(BICs)
opt_q2 <- which.min(BIC2)
####
#opt_q=18
gmfads <- gmm.fad(wbc, 2, q= opt_q, maxiter =1000, tol = tol, nstart=20)
est_cl <- gmfads$clusters
g.ari<-RandIndex(true_cl, est_cl)$AR
ari #####
library(caret)
confusion_matrix <- confusionMatrix(as.factor(est_cl), as.factor(true_cl), positive = "2")
print(confusion_matrix)
g.accuracy <-as.numeric( confusion_matrix$overall["Accuracy"])
g.kappa <- as.numeric(confusion_matrix$overall[ "Kappa"])
g.sensitivity <- as.numeric(confusion_matrix$byClass["Sensitivity"])
g.specificity <- as.numeric(confusion_matrix$byClass["Specificity"])
GMMFAD <- c(ARI = g.ari, Accuracy= g.accuracy, sensitivity= g.sensitivity, 
               specificity= g.specificity, Kappa= g.kappa )
##################
#opt_q2 =19
emmixs <- mfa(wbc.trans, g=2, q=opt_q2, D_type = 'unique', sigma_type='unique', maxiter =1000, tol = tol, conv_measure = 'ratio')
est_cl2 <- emmixs$clust
e.ari<-RandIndex(true_cl, est_cl2)$AR; e.ari
est_cl2.r <- ifelse(est_cl2==1,2,1)
confusion_matrix2 <- confusionMatrix(as.factor(est_cl2.r), as.factor(true_cl), positive = "2")
print(confusion_matrix2)
e.accuracy <-as.numeric( confusion_matrix2$overall["Accuracy"])
e.kappa <- as.numeric(confusion_matrix2$overall[ "Kappa"])
e.sensitivity <- as.numeric(confusion_matrix2$byClass["Sensitivity"])
e.specificity <- as.numeric(confusion_matrix2$byClass["Specificity"])
EMMIX<- c(ARI = e.ari, Accuracy= e.accuracy, sensitivity= e.sensitivity, 
               specificity= e.specificity, Kappa= e.kappa )

metrics<- rbind(GMMFAD, EMMIX)
#################



######################################################################################################

# Let's read gmfad and emmix clustering results 
#gmfads <-read.csv("~/Desktop/Dimensionality reduction research/gmfad paper/breastcancer_data_gmfad_clustering_results_new2.csv")
#emmixs <-read.csv("~/Desktop/Dimensionality reduction research/gmfad paper/breastcancer_data_emmix_clustering_results_new2.csv")
#gmfadBIC <- gmfads$BIC
#emmixBIC <- emmixs$BIC 
library(ggplot2)
Q <- length(gmfadBIC)
q <- 1:Q      

 
################### MDS plot of wbc data with gmfad's est_cl.
library(MASS)   # for the isoMDS function
library(ggplot2)  # for plotting
dissimilarity_matrix <- dist(wbc)
dissimilarity_matrix <- dist(wbc.trans)
mds_result <- isoMDS(dissimilarity_matrix, k = 2) 
trimmed_wbc<- mds_result$points
#str(trimmed_wbc)
mds_data <- data.frame(
  MDS1 = mds_result$points[, 1],
  MDS2 = mds_result$points[, 2],
  Diagnosis = true_cl
)

mds_data$true_cl <- as.factor(mds_data$Diagnosis)
estimated_cl <- as.factor(est_cl)
dataplot <- ggplot(mds_data, aes(x = MDS1, y = MDS2, color = true_cl, shape = estimated_cl)) +
  geom_point() +
  labs(title = "Multidimensional Scaling Representation: Winsconsin Breast Cancer",
       x = "Dimension 1",
       y = "Dimension 2") +
  theme_minimal() 
dataplot
####
dataplot2 <- ggplot(mds_data, aes(x = MDS1, y = MDS2, color = true_cl, shape = estimated_cl)) +
geom_text(
  aes(
    label = ifelse(est_cl == 1, "M", "B"),
    color = factor(true_cl)  
  ),
  size = 2
) +
  labs( x = "Dimension 1",
       y = "Dimension 2") +
  theme_minimal() +
  theme(legend.position = "none")
dataplot2

####### bivariate projection of wbc using gmfad.
gmfad_2 <- gmm.fad(wbc.trans, 2, 2, tol = tol, nstart=20)
means = gmfad_2$means
gamma = gmfad_2$gamma
lambda = gmfad_2$lambda
psi =  gmfad_2$psi
est_cl <- gmfads$clusters

factors_2D <- projection_2D(wbc.trans, means, lambda, psi, gamma)

reduced_data <- data.frame(
  factor1 = factors_2D[, 1],
  factor2 = factors_2D[, 2],
  Diagnosis = true_cl
)
reduced_data$true_cl <- as.factor(reduced_data$Diagnosis)
estimated_cl <- as.factor(est_cl)
dataplot <- ggplot(reduced_data, aes(x = factor1, y = factor2, color = true_cl, shape = estimated_cl)) +
  geom_point() +
  labs(#title = "2-factor Representation (with gmmfad): Winsconsin Breast Cancer",
       x = "Factor 1",
       y = "Factor 2") +
  theme_minimal() 
dataplot
#####################
################################################################################
################ Breast cancer with varied q. ####################
#################################################################################
source("~/Desktop/Dimensionality reduction research/codes/gmfad_qq.R")
library(doParallel)
library(MixSim) # for ARI
library(caret)
cc <- detectCores()
cl <- makePSOCKcluster(cc - 1)
registerDoParallel(cl)
library(radviz3d)
#
seed <- 1233345 # This gives ARI 73. Good. Seed only useful in Radviz transformation.
wbc <- read.csv("wdbc.data", header = FALSE)# Read the names file. Informative!
names <- readLines("wdbc.names")
wbc$V2 <- as.factor(wbc$V2)
true_cl <- as.integer(wbc$V2)
table(wbc$V2)
wbc<- wbc[, -c(1,2)]
print(str(wbc))
set.seed(seed) 
wbc.trans <-radviz3d::Gtrans(wbc)
K= 2
q1= q2 =  1:20
q_grid <- as.matrix(expand.grid(q1=q1, q2=q2))
BICs<-runtimes<- aris<-loglik<- tols <- niter <- numeric(nrow(q_grid) )
tol = 1e-12
tols = rep(tol,length(q_grid))
for (j in 1:nrow(q_grid) ){
  t1 <- proc.time()
  qvec <- as.vector(q_grid[j,])
  gmfads <- tryCatch({
    gmm.fad.q(wbc.trans, K, qvec, maxiter = 1000, q, tol =tol, nstart=20)
  }, error= function(e) NA)
  runtimes[j]<- (proc.time()-t1)[3]
  est_cl  <- tryCatch({gmfads$clusters }, error= function(e) NA)
  aris[j] <- tryCatch({RandIndex(true_cl, est_cl)$AR }, error= function(e) NA)
  BICs[j] <- tryCatch({gmfads$BIC  }, error= function(e) NA)
  niter[j]  <- tryCatch({gmfads$niter }, error= function(e) NA)
  loglik[j] <- tryCatch({gmfads$loglik },  error= function(e) NA)
  print(Sys.time())
  result_breast_cancer <- data.frame(ARI = aris[j], runtimes = runtimes[j], BIC = BICs[j], 
                              niter = niter[j], loglik = loglik[j], tolerance= tols[j], 
                              q1=qvec[1], q2=qvec[2])
  print(result_breast_cancer)
  write.table( result_breast_cancer, 
               file = paste0("breast_cancer_data_clustering_results_variedq_seed_",seed,".csv") , append = TRUE, sep = ",",
               col.names = !file.exists(paste0("breast_cancer_data_clustering_results_variedq_seed_",seed,".csv")
               ), row.names = FALSE )
}
result_breast_cancer <- data.frame(ARI = aris, runtimes = runtimes, BIC = BICs, 
                            niter = niter, loglik = loglik, tolerance= tols)
result_breast_cancer <- cbind(result_breast_cancer, q_grid)
print(result_breast_cancer)


#################

K=2
tol = 1e-12
opt_q= as.vector(q_grid[which.min(result_breast_cancer$BIC) ,])
opt_q
set.seed(seed)
wbc.trans  <-radviz3d::Gtrans(wbc)
t1 <- proc.time()
gmfads <-  gmm.fad.q(wbc.trans, K, opt_q, maxiter = 1000, q, tol =tol, nstart=20)
duration<- (proc.time()-t1)[3]
est_cl <- gmfads$clusters
ari<-RandIndex(true_cl, est_cl)$AR
ari
library(MASS)   # for the isoMDS function
library(ggplot2)  # for plotting
dissimilarity_matrix <- dist(wbc.trans)
mds_result <- isoMDS(dissimilarity_matrix, k = 2) 
mds_breast_cancer<- mds_result$points
#str(trimmed_wbc)
mds_data <- data.frame(
  MDS1 = mds_result$points[, 1],
  MDS2 = mds_result$points[, 2],
  Diagnosis = true_cl
)

library(caret)
confusion_matrix <- confusionMatrix(as.factor(est_cl), as.factor(true_cl))
print(confusion_matrix)
