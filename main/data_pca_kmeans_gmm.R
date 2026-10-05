## ============================================================================
## PCA-based clustering baselines with data-driven PC selection
##
## Baselines:
##   1. PCA + k-means
##      PC number selected by maximum average silhouette width
##
##   2. PCA + GMM
##      PC number selected by V-fold cross-validated validation log-likelihood
##
## Outputs:
##   - selected number of PCs
##   - ARI
##   - elapsed runtime in seconds
##   - completion/convergence status
##
## Important:
##   - True class labels are NOT used to select the number of PCs.
##   - ARI is used only for external evaluation after tuning.
##   - Runtime includes PC selection and final model fitting.
## ============================================================================


## ---------------------------------------------------------------------------
## 0. Packages
## ---------------------------------------------------------------------------

required_pkgs <- c(
  "MixSim",
  "mclust",
  "spls",
  "radviz3d",
  "cluster"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg)
  }
}

library(MixSim)
library(mclust)
library(spls)
library(cluster)


## ---------------------------------------------------------------------------
## 1. ARI
## ---------------------------------------------------------------------------

compute_ari <- function(true_cl, pred_cl) {
  
  ri <- try(
    MixSim::RandIndex(true_cl, pred_cl),
    silent = TRUE
  )
  
  if (!inherits(ri, "try-error")) {
    
    if (is.list(ri) && "AR" %in% names(ri)) {
      return(as.numeric(ri$AR))
    }
    
    if (is.list(ri)) {
      
      nm <- names(ri)
      
      id <- grep(
        "adj|adjust|AR",
        nm,
        ignore.case = TRUE
      )
      
      if (length(id) > 0) {
        return(as.numeric(ri[[id[1]]]))
      }
    }
  }
  
  mclust::adjustedRandIndex(
    true_cl,
    pred_cl
  )
}


## ---------------------------------------------------------------------------
## 2. Timing helper
## ---------------------------------------------------------------------------

measure_fun <- function(fun) {
  
  value <- NULL
  error_message <- NA_character_
  status <- "completed"
  
  t0 <- proc.time()
  
  value <- tryCatch(
    fun(),
    error = function(e) {
      
      status <<- "failed"
      
      error_message <<-
        conditionMessage(e)
      
      NULL
    }
  )
  
  runtime_sec <-
    as.numeric(
      (proc.time() - t0)["elapsed"]
    )
  
  list(
    value = value,
    runtime_sec = runtime_sec,
    status = status,
    error_message = error_message
  )
}


## ---------------------------------------------------------------------------
## 3. Basic data preparation
## ---------------------------------------------------------------------------

prepare_X <- function(X) {
  
  X <- as.matrix(X)
  
  storage.mode(X) <-
    "double"
  
  ## Remove zero-variance variables.
  sds <-
    apply(
      X,
      2,
      sd,
      na.rm = TRUE
    )
  
  keep <-
    is.finite(sds) &
    sds > 0
  
  X <-
    X[, keep, drop = FALSE]
  
  if (ncol(X) < 1) {
    stop(
      "No non-zero-variance features remain."
    )
  }
  
  X
}


## ---------------------------------------------------------------------------
## 4. PCA on complete dataset
## ---------------------------------------------------------------------------

make_pca_scores <- function(
    X,
    max_pcs) {
  
  X <- prepare_X(X)
  
  max_pcs_eff <-
    min(
      max_pcs,
      nrow(X) - 1,
      ncol(X)
    )
  
  if (max_pcs_eff < 1) {
    stop(
      "No principal components are available."
    )
  }
  
  pca <-
    prcomp(
      X,
      center = TRUE,
      scale. = TRUE,
      rank. = max_pcs_eff
    )
  
  Z <-
    pca$x[
      ,
      seq_len(max_pcs_eff),
      drop = FALSE
    ]
  
  list(
    scores = Z,
    pca = pca,
    max_pcs = max_pcs_eff,
    n_features_used = ncol(X)
  )
}


## ---------------------------------------------------------------------------
## 5. PCA projection for cross-validation
## ---------------------------------------------------------------------------
## PCA is fitted using training observations only.
## Validation observations are standardized and projected using
## training means, SDs, and PCA rotation.

fit_project_pca <- function(
    X_train,
    X_valid,
    max_pcs) {
  
  X_train <- as.matrix(X_train)
  X_valid <- as.matrix(X_valid)
  
  storage.mode(X_train) <- "double"
  storage.mode(X_valid) <- "double"
  
  train_means <-
    colMeans(X_train)
  
  train_sds <-
    apply(
      X_train,
      2,
      sd
    )
  
  keep <-
    is.finite(train_sds) &
    train_sds > 0
  
  X_train <-
    X_train[, keep, drop = FALSE]
  
  X_valid <-
    X_valid[, keep, drop = FALSE]
  
  train_means <-
    train_means[keep]
  
  train_sds <-
    train_sds[keep]
  
  X_train_scaled <-
    sweep(
      X_train,
      2,
      train_means,
      "-"
    )
  
  X_train_scaled <-
    sweep(
      X_train_scaled,
      2,
      train_sds,
      "/"
    )
  
  X_valid_scaled <-
    sweep(
      X_valid,
      2,
      train_means,
      "-"
    )
  
  X_valid_scaled <-
    sweep(
      X_valid_scaled,
      2,
      train_sds,
      "/"
    )
  
  max_pcs_eff <-
    min(
      max_pcs,
      nrow(X_train_scaled) - 1,
      ncol(X_train_scaled)
    )
  
  pca <-
    prcomp(
      X_train_scaled,
      center = FALSE,
      scale. = FALSE,
      rank. = max_pcs_eff
    )
  
  rotation <-
    pca$rotation[
      ,
      seq_len(max_pcs_eff),
      drop = FALSE
    ]
  
  train_scores <-
    X_train_scaled %*% rotation
  
  valid_scores <-
    X_valid_scaled %*% rotation
  
  list(
    train_scores = train_scores,
    valid_scores = valid_scores,
    max_pcs = max_pcs_eff
  )
}


## ===========================================================================
## 6. PCA + k-means
##    PC selection by average silhouette width
## ===========================================================================

run_pca_kmeans <- function(
    X,
    true_cl,
    dataset_name,
    pc_grid = 2:20,
    nstart = 20,
    seed = 1234456) {
  
  true_cl <-
    as.factor(true_cl)
  
  X <-
    as.matrix(X)
  
  keep_rows <-
    complete.cases(X) &
    !is.na(true_cl)
  
  X <-
    X[
      keep_rows,
      ,
      drop = FALSE
    ]
  
  true_cl <-
    droplevels(
      true_cl[keep_rows]
    )
  
  K <-
    nlevels(true_cl)
  
  measured <-
    measure_fun(
      function() {
        
        ## PCA only once up to maximum candidate dimension.
        pca_out <-
          make_pca_scores(
            X,
            max_pcs = max(pc_grid)
          )
        
        Z <-
          pca_out$scores
        
        valid_grid <-
          pc_grid[
            pc_grid <= ncol(Z)
          ]
        
        if (length(valid_grid) == 0) {
          stop(
            "No valid PC dimensions."
          )
        }
        
        silhouette_scores <-
          rep(
            NA_real_,
            length(valid_grid)
          )
        
        ## ---------------------------------------------------------------
        ## Select PC number
        ## ---------------------------------------------------------------
        
        for (i in seq_along(valid_grid)) {
          
          m <-
            valid_grid[i]
          
          Zm <-
            Z[
              ,
              seq_len(m),
              drop = FALSE
            ]
          
          set.seed(
            seed + m
          )
          
          km <-
            kmeans(
              Zm,
              centers = K,
              nstart = nstart,
              iter.max = 1000
            )
          
          sil <-
            cluster::silhouette(
              km$cluster,
              dist(Zm)
            )
          
          silhouette_scores[i] <-
            mean(
              sil[, "sil_width"]
            )
        }
        
        best_id <-
          which.max(
            silhouette_scores
          )
        
        selected_pcs <-
          valid_grid[best_id]
        
        ## ---------------------------------------------------------------
        ## Final fit
        ## ---------------------------------------------------------------
        
        Z_final <-
          Z[
            ,
            seq_len(selected_pcs),
            drop = FALSE
          ]
        
        set.seed(seed)
        
        final_fit <-
          kmeans(
            Z_final,
            centers = K,
            nstart = nstart,
            iter.max = 1000
          )
        
        list(
          fit = final_fit,
          selected_pcs = selected_pcs,
          selection_score =
            silhouette_scores[best_id],
          pc_grid = valid_grid,
          silhouette_scores =
            silhouette_scores,
          n_features_used =
            pca_out$n_features_used
        )
      }
    )
  
  
  if (
    measured$status == "failed" ||
    is.null(measured$value)
  ) {
    
    return(
      data.frame(
        dataset = dataset_name,
        method = "PCA+k-means",
        K = K,
        selected_pcs = NA_integer_,
        selection_criterion =
          "Average silhouette width",
        selection_score = NA_real_,
        ARI = NA_real_,
        runtime_sec =
          measured$runtime_sec,
        status = "failed",
        error_message =
          measured$error_message,
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  fit <-
    measured$value$fit
  
  ari <-
    compute_ari(
      true_cl,
      fit$cluster
    )
  
  
  ## Hartigan-Wong diagnostic.
  if (
    !is.null(fit$ifault) &&
    fit$ifault != 0
  ) {
    
    fit_status <-
      paste0(
        "completed (ifault=",
        fit$ifault,
        ")"
      )
    
  } else {
    
    fit_status <-
      "completed"
  }
  
  
  data.frame(
    dataset = dataset_name,
    method = "PCA+k-means",
    K = K,
    selected_pcs =
      measured$value$selected_pcs,
    selection_criterion =
      "Average silhouette width",
    selection_score =
      measured$value$selection_score,
    ARI = ari,
    runtime_sec =
      measured$runtime_sec,
    status = fit_status,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}


## ===========================================================================
## 7. PCA + GMM
##    PC selection by cross-validated validation log-likelihood
## ===========================================================================

run_pca_gmm <- function(
    X,
    true_cl,
    dataset_name,
    pc_grid = 2:20,
    nfolds = 5,
    seed = 1234456) {
  
  true_cl <-
    as.factor(true_cl)
  
  X <-
    as.matrix(X)
  
  keep_rows <-
    complete.cases(X) &
    !is.na(true_cl)
  
  X <-
    X[
      keep_rows,
      ,
      drop = FALSE
    ]
  
  true_cl <-
    droplevels(
      true_cl[keep_rows]
    )
  
  K <-
    nlevels(true_cl)
  
  
  measured <-
    measure_fun(
      function() {
        
        n <-
          nrow(X)
        
        ## ---------------------------------------------------------------
        ## Construct folds independently of class labels
        ## ---------------------------------------------------------------
        
        set.seed(seed)
        
        fold_id <-
          sample(
            rep(
              seq_len(nfolds),
              length.out = n
            )
          )
        
        
        cv_loglik <-
          rep(
            NA_real_,
            length(pc_grid)
          )
        
        
        ## ---------------------------------------------------------------
        ## Evaluate each PC dimension
        ## ---------------------------------------------------------------
        
        for (j in seq_along(pc_grid)) {
          
          m <-
            pc_grid[j]
          
          total_loglik <-
            0
          
          total_n <-
            0
          
          failed <-
            FALSE
          
          
          for (fold in seq_len(nfolds)) {
            
            train_id <-
              fold_id != fold
            
            valid_id <-
              fold_id == fold
            
            
            pca_cv <-
              try(
                fit_project_pca(
                  X_train =
                    X[train_id, , drop = FALSE],
                  X_valid =
                    X[valid_id, , drop = FALSE],
                  max_pcs = m
                ),
                silent = TRUE
              )
            
            
            if (
              inherits(
                pca_cv,
                "try-error"
              ) ||
              pca_cv$max_pcs < m
            ) {
              
              failed <- TRUE
              
              break
            }
            
            
            Z_train <-
              pca_cv$train_scores[
                ,
                seq_len(m),
                drop = FALSE
              ]
            
            Z_valid <-
              pca_cv$valid_scores[
                ,
                seq_len(m),
                drop = FALSE
              ]
            
            
            gmm_fit <-
              try(
                mclust::Mclust(
                  Z_train,
                  G = K,
                  verbose = FALSE
                ),
                silent = TRUE
              )
            
            
            if (
              inherits(
                gmm_fit,
                "try-error"
              ) ||
              is.null(
                gmm_fit$parameters
              )
            ) {
              
              failed <- TRUE
              
              break
            }
            
            
            logdens <-
              try(
                mclust::dens(
                  modelName =
                    gmm_fit$modelName,
                  data =
                    Z_valid,
                  parameters =
                    gmm_fit$parameters,
                  logarithm = TRUE
                ),
                silent = TRUE
              )
            
            
            if (
              inherits(
                logdens,
                "try-error"
              ) ||
              any(
                !is.finite(logdens)
              )
            ) {
              
              failed <- TRUE
              
              break
            }
            
            
            total_loglik <-
              total_loglik +
              sum(logdens)
            
            total_n <-
              total_n +
              length(logdens)
          }
          
          
          if (
            !failed &&
            total_n > 0
          ) {
            
            ## Mean validation log-likelihood
            ## per observation.
            cv_loglik[j] <-
              total_loglik /
              total_n
          }
        }
        
        
        valid <-
          is.finite(
            cv_loglik
          )
        
        
        if (!any(valid)) {
          
          stop(
            "All PCA+GMM cross-validation fits failed."
          )
        }
        
        
        best_id <-
          which.max(
            cv_loglik
          )
        
        selected_pcs <-
          pc_grid[best_id]
        
        
        ## ---------------------------------------------------------------
        ## Final PCA on full dataset
        ## ---------------------------------------------------------------
        
        pca_final <-
          make_pca_scores(
            X,
            max_pcs = selected_pcs
          )
        
        
        Z_final <-
          pca_final$scores[
            ,
            seq_len(selected_pcs),
            drop = FALSE
          ]
        
        
        ## ---------------------------------------------------------------
        ## Final GMM
        ## ---------------------------------------------------------------
        
        final_fit <-
          mclust::Mclust(
            Z_final,
            G = K,
            verbose = FALSE
          )
        
        
        list(
          fit = final_fit,
          selected_pcs =
            selected_pcs,
          selection_score =
            cv_loglik[best_id],
          pc_grid =
            pc_grid,
          cv_loglik =
            cv_loglik
        )
      }
    )
  
  
  if (
    measured$status == "failed" ||
    is.null(measured$value)
  ) {
    
    return(
      data.frame(
        dataset = dataset_name,
        method = "PCA+GMM",
        K = K,
        selected_pcs = NA_integer_,
        selection_criterion =
          "CV log-likelihood",
        selection_score = NA_real_,
        ARI = NA_real_,
        runtime_sec =
          measured$runtime_sec,
        status = "failed",
        error_message =
          measured$error_message,
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  fit <-
    measured$value$fit
  
  
  if (
    is.null(
      fit$classification
    )
  ) {
    
    ari <-
      NA_real_
    
    fit_status <-
      "failed"
    
  } else {
    
    ari <-
      compute_ari(
        true_cl,
        fit$classification
      )
    
    fit_status <-
      "completed"
  }
  
  
  data.frame(
    dataset = dataset_name,
    method = "PCA+GMM",
    K = K,
    selected_pcs =
      measured$value$selected_pcs,
    selection_criterion =
      "CV log-likelihood",
    selection_score =
      measured$value$selection_score,
    ARI = ari,
    runtime_sec =
      measured$runtime_sec,
    status = fit_status,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}


## ===========================================================================
## 8. Wrapper for both baselines
## ===========================================================================

run_pca_baselines <- function(
    X,
    true_cl,
    dataset_name,
    pc_grid = 2:20,
    nstart = 20,
    nfolds = 5,
    seed = 1234456) {
  
  res_kmeans <-
    run_pca_kmeans(
      X = X,
      true_cl = true_cl,
      dataset_name = dataset_name,
      pc_grid = pc_grid,
      nstart = nstart,
      seed = seed
    )
  
  res_gmm <-
    run_pca_gmm(
      X = X,
      true_cl = true_cl,
      dataset_name = dataset_name,
      pc_grid = pc_grid,
      nfolds = nfolds,
      seed = seed
    )
  
  rbind(
    res_kmeans,
    res_gmm
  )
}


## ===========================================================================
## 9. Data applications
## ===========================================================================

setwd("~/Desktop/gmfad paper")


## ---------------------------------------------------------------------------
## Wisconsin breast cancer
## ---------------------------------------------------------------------------

wbc <-
  read.csv(
    "wdbc.data",
    header = FALSE
  )

wbc$V2 <-
  as.factor(
    wbc$V2
  )

true_cl_wbc <-
  wbc$V2

wbc_x <-
  wbc[, -c(1, 2)]

wbc_trans <-
  radviz3d::Gtrans(
    wbc_x
  )


res_wbc <-
  run_pca_baselines(
    X = wbc_trans,
    true_cl = true_cl_wbc,
    dataset_name =
      "Wisconsin breast cancer",
    pc_grid = 2:20,
    nstart = 20,
    nfolds = 5,
    seed = 1234456
  )


## ---------------------------------------------------------------------------
## Lymphoma
## ---------------------------------------------------------------------------

data(
  lymphoma,
  package = "spls"
)

lymph_x <-
  lymphoma$x

true_cl_lymph <-
  as.factor(
    lymphoma$y + 1
  )


res_lymph <-
  run_pca_baselines(
    X = lymph_x,
    true_cl = true_cl_lymph,
    dataset_name =
      "Lymphoma gene expression",
    pc_grid = 2:20,
    nstart = 20,
    nfolds = 5,
    seed = 1234456
  )


## ---------------------------------------------------------------------------
## Human brain scRNA-seq
## ---------------------------------------------------------------------------

load(
  "darmanis-adult-neuron-astrocyte-200HVG-GDT.rda"
)

scrna_x <-
  as.matrix(
    darmanis_adult_2class_df[, -1]
  )

true_cl_scrna <-
  as.factor(
    darmanis_adult_2class_df$cell_type
  )


res_scrna <-
  run_pca_baselines(
    X = scrna_x,
    true_cl = true_cl_scrna,
    dataset_name =
      "Human brain scRNA-seq",
    pc_grid = 2:20,
    nstart = 20,
    nfolds = 5,
    seed = 1234456
  )


## ===========================================================================
## 10. Combine results
## ===========================================================================

baseline_results <-
  rbind(
    res_wbc,
    res_lymph,
    res_scrna
  )

print(
  baseline_results
)


## ===========================================================================
## 11. Rounded output for manuscript
## ===========================================================================

baseline_results_print <-
  baseline_results

baseline_results_print$ARI <-
  round(
    baseline_results_print$ARI,
    3
  )

baseline_results_print$selection_score <-
  round(
    baseline_results_print$selection_score,
    4
  )

baseline_results_print$runtime_sec <-
  round(
    baseline_results_print$runtime_sec,
    2
  )

print(
  baseline_results_print
)


## ===========================================================================
## 12. Save
## ===========================================================================

write.csv(
  baseline_results,
  "real_data_pca_kmeans_gmm_tuned_results.csv",
  row.names = FALSE
)
