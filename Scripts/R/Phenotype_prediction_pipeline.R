library(data.table)
library(caret)
library(ranger)
library(glmnet)
library(parallel)
library(doParallel)

# --- SÉLECTEUR DE VARIABLES  ---
get_snp_selection <- function(type, all_scores_dt, ph, m, ctx, topK, X_names) {
  if (type == "random") return(sample(X_names, topK))

  if (type == "score") {
    snps <- all_scores_dt[phenotype == ph & method == m & context == ctx][order(rank)][1:topK, SNP]
    return(snps)
  }

  if (type == "hybrid") {
    s1 <- all_scores_dt[phenotype == ph & method == "FarmCPU"][order(rank)][1:floor(topK/2), SNP]
    s2 <- all_scores_dt[phenotype == ph & method == "RF"][order(rank)][1:floor(topK/2), SNP]
    return(unique(c(s1, s2)))
  }
  return(character(0))
}

# --- MOTEUR CENTRAL ---
core_nested_cv <- function(X, y, model_type, model_params, K_outer=5, K_inner=5, seed=123) {

  set.seed(seed)
  outer_folds <- createFolds(y, k = K_outer, returnTrain = TRUE)
  results_dt <- data.table()

  for (f in seq_len(K_outer)) {

    tr_idx <- outer_folds[[f]]
    te_idx <- setdiff(seq_along(y), tr_idx)

    X_tr <- X[tr_idx, , drop = FALSE]
    y_tr <- y[tr_idx]
    X_te <- X[te_idx, , drop = FALSE]
    y_te <- y[te_idx]

    y_pred <- NULL
    best_params <- list(mtry = NA, node = NA, lambda = NA)

    # --- MODEL RANDOM FOREST ---
    if (model_type == "RF") {
      # 1. Tuning (CV Interne)
      inner_folds <- createFolds(y_tr, k = K_inner, returnTrain = TRUE)
      grid_perf <- as.data.table(copy(model_params$rf_grid))
      grid_perf[, RMSE_cv := NA_real_]

      for (g in seq_len(nrow(grid_perf))) {
        rmse_in <- numeric(K_inner)
        for (k_in in seq_len(K_inner)) {
          in_tr <- inner_folds[[k_in]]
          in_te <- setdiff(seq_along(y_tr), in_tr)

          mod_in <- ranger(
            x = as.data.frame(X_tr[in_tr, ]), y = y_tr[in_tr],
            num.trees = grid_perf$num.trees[g],
            mtry = min(grid_perf$mtry[g], ncol(X_tr)),
            min.node.size = grid_perf$min.node.size[g],
            num.threads = 1, verbose = FALSE
          )
          p_in <- predict(mod_in, as.data.frame(X_tr[in_te, ]))$predictions
          rmse_in[k_in] <- sqrt(mean((y_tr[in_te] - p_in)^2))
        }
        grid_perf$RMSE_cv[g] <- mean(rmse_in)
      }

      # Fit Final
      best <- grid_perf[which.min(RMSE_cv)]
      best_params$mtry <- best$mtry
      best_params$node <- best$min.node.size

      final_mod <- ranger(
        x = as.data.frame(X_tr), y = y_tr,
        num.trees = best$num.trees,
        mtry = min(best$mtry, ncol(X_tr)),
        min.node.size = best$min.node.size,
        importance = "permutation", num.threads = 1, verbose = FALSE
      )
      y_pred <- predict(final_mod, as.data.frame(X_te))$predictions
    }

    # --- MODEL RIDGE ---
    else if (model_type == "ridge") {
      cv_fit <- cv.glmnet(X_tr, y_tr, alpha = 0, nfolds = K_inner, type.measure = "mse")
      best_params$lambda <- cv_fit$lambda.min
      y_pred <- as.vector(predict(cv_fit, newx = X_te, s = "lambda.min"))
    }

    # --- CALCUL DES MÉTRIQUES (R2 & RMSE) ---
    sse <- sum((y_te - y_pred)^2)
    sst <- sum((y_te - mean(y_te))^2)

    r2_val <- 1 - (sse / sst)
    rmse_val <- sqrt(mean((y_te - y_pred)^2))


    res_row <- data.table(
      outer_fold = f,
      R2 = r2_val,
      RMSE = rmse_val,
      best_mtry = best_params$mtry,
      best_node = best_params$node,
      best_lambda = best_params$lambda
    )

    results_dt <- rbind(results_dt, res_row)
  }

  return(results_dt)
}



run_parallel_pipeline <- function(
    topK_vector, all_scores_dt, phenotypes, methods, contexts,
    X_full, y_full_list,
    model_type = "RF",
    rf_grid = NULL,
    n_cv_repeats = 5,
    n_random_reps = 3,
    K_outer = 5, K_inner = 5,
    n_cores = 4,
    seed_start = 123,
    log_file = "pipeline_progress.txt"
) {

  # --- 1. PRÉPARATION DES TÂCHES ---
  tasks <- list()
  X_cols <- colnames(X_full)
  for (K in topK_vector) {
    for (ph in phenotypes) {
      for (rr in seq_len(n_random_reps)) {
        for (cvr in seq_len(n_cv_repeats)) {
          tasks[[length(tasks) + 1]] <- list(task_id = "random", K = K, pheno = ph,
                                             method = "Random", context = "Baseline",
                                             seed = seed_start + (cvr * 1000) + rr, cv_rep = cvr)
        }
      }
      for (m in methods) {
        for (ctx in contexts) {
          if(nrow(all_scores_dt[phenotype==ph & method==m & context==ctx]) > 0) {
            for (cvr in seq_len(n_cv_repeats)) {
              tasks[[length(tasks) + 1]] <- list(task_id = "score", K = K, pheno = ph,
                                                 method = m, context = ctx,
                                                 seed = seed_start + (cvr * 1000), cv_rep = cvr)
            }
          }
        }
      }
    }
  }

  n_tasks <- length(tasks)
  cat(paste0(Sys.time(), " | Début : ", n_tasks, " tâches prévues.\n"), file = log_file)

  # ---  EXÉCUTION PARALLÈLE ---
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  final_results <- foreach(
    i = seq_along(tasks),
    .packages = c("data.table", "ranger", "caret", "glmnet"),
    .export = c("core_nested_cv", "get_snp_selection"),
    .combine = rbind,
    .errorhandling = "remove"
  ) %dopar% {

    task <- tasks[[i]]

    if (i %% 5 == 0) {
      cat(paste0(Sys.time(), " | Progress: ", round(i/n_tasks*100, 1), "% (Task ", i, "/", n_tasks, ")\n"),
          file = log_file, append = TRUE)
    }

    sel_snps <- get_snp_selection(
      type = ifelse(task$task_id == "random", "random", "score"),
      all_scores_dt = all_scores_dt,
      ph = task$pheno, m = task$method, ctx = task$context,
      topK = task$K, X_names = X_cols
    )

    sel_snps <- sel_snps[sel_snps %in% X_cols]
    if (length(sel_snps) < 2) return(NULL)

    perf <- core_nested_cv(
      X = X_full[, sel_snps, drop = FALSE],
      y = y_full_list[[task$pheno]],
      model_type = model_type,
      model_params = list(rf_grid = rf_grid),
      K_outer = K_outer, K_inner = K_inner,
      seed = task$seed
    )

    perf[, `:=`(phenotype = task$pheno, method = task$method, context = task$context,
                topK = task$K, cv_repeat = task$cv_rep, n_snps = length(sel_snps))]

    return(perf)
  }

  stopCluster(cl)
  cat(paste0(Sys.time(), " | Terminé !\n"), file = log_file, append = TRUE)

  if(!is.null(final_results)) {
    final_results[, topK_factor := factor(topK, levels = sort(unique(topK)))]
  }

  return(final_results)
}
