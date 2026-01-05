library(pls)
library(parallel)
library(doParallel)
library(data.table)
library(glmnet)

run_ridge_nested_cv <- function(
    X,
    y,
    K_outer = 5,
    K_inner = 5,
    seed = 123,
    standardize = TRUE
) {
  stopifnot(nrow(X) == length(y))

  set.seed(seed)

  library(caret)
  library(glmnet)

  # Création des folds externes
  outer_folds <- createFolds(y, k = K_outer, returnTrain = TRUE)

  results <- data.frame(
    outer_fold = integer(),
    lambda_opt = numeric(),
    R2 = numeric(),
    RMSE = numeric()
  )

  for (f in seq_len(K_outer)) {

    # Split externe
    train_idx <- outer_folds[[f]]
    test_idx  <- setdiff(seq_along(y), train_idx)

    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx, , drop = FALSE]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]

    # CV interne pour le choix du lambda
    cv_inner <- cv.glmnet(
      x = X_train,
      y = y_train,
      alpha = 0,                # Ridge
      nfolds = K_inner,
      standardize = standardize
    )

    lambda_opt <- cv_inner$lambda.min

    # Modèle final entraîné sur tout le train externe
    ridge_final <- glmnet(
      x = X_train,
      y = y_train,
      alpha = 0,
      lambda = lambda_opt,
      standardize = standardize
    )

    # Prédictions sur le test externe
    y_pred <- as.vector(predict(ridge_final, newx = X_test))

    #R2_test <- cor(y_test, y_pred)^2
    R2_test <- 1 - mean((y_test - y_pred)^2) / var(y_test)
    RMSE_test <- sqrt(mean((y_test - y_pred)^2))

    results <- rbind(
      results,
      data.frame(
        outer_fold = f,
        lambda_opt = lambda_opt,
        R2 = R2_test,
        RMSE = RMSE_test
      )
    )
  }

  return(results)
}

################################################################################
run_ridge_from_scores <- function(
    all_scores_dt,
    phenotypes,
    methods,
    contexts,
    topK,
    X_full,
    y_full,
    K_outer = 5,
    K_inner = 5,
    seed = 123,
    standardize = TRUE
) {

  library(data.table)
  cat(">>> ENTERING run_model_from_scores <<<\n")

  set.seed(seed)   # pour reproductibilité
  res_all <- list()

  for (ph in phenotypes) {
    ## --- Ajouter le témoin aléatoire ---
    snp_random <- sample(colnames(X_full), topK)
    X_rand <- X_full[, snp_random, drop = FALSE]

    perf_rand <- run_ridge_nested_cv(
      X = X_rand,
      y = y_full,
      K_outer = K_outer,
      K_inner = K_inner,
      seed = seed + 1,  # décaler le seed pour varier
      standardize = standardize
    )

    perf_rand_dt <- as.data.table(perf_rand)
    perf_rand_dt[, `:=`(
      phenotype = ph,
      method = "Random",
      context = "temoin",
      topK = topK,
      n_snps = ncol(X_rand)
    )]

    res_all[[length(res_all) + 1]] <- perf_rand_dt

    for (m in methods) {
      for (ctx in contexts) {

        message(
          "Running ridge | phenotype = ", ph,
          " | method = ", m,
          " | context = ", ctx
        )

        cat(
          "Running ridge | phenotype = ", ph,
          " | method = ", m,
          " | context = ", ctx
        )

        ## Sélection des TOP K SNPs
        snp_sel <- all_scores_dt[
          phenotype == ph &
            method == m &
            context == ctx
        ][order(rank)][1:topK, SNP]

        if (length(snp_sel) < 2) next

        ## Sous-matrice X
        X_sub <- X_full[, colnames(X_full) %in% snp_sel, drop = FALSE]

        ## Run ridge nested CV
        perf <- run_ridge_nested_cv(
          X = X_sub,
          y = y_full,
          K_outer = K_outer,
          K_inner = K_inner,
          seed = seed,
          standardize = standardize
        )

        perf_dt <- as.data.table(perf)
        perf_dt[, `:=`(
          phenotype = ph,
          method = m,
          context = ctx,
          topK = topK,
          n_snps = ncol(X_sub)
        )]

        res_all[[length(res_all) + 1]] <- perf_dt
      }
    }
  }

  return(rbindlist(res_all, fill = TRUE))
}



################################################################################

run_ridge_model <- function(
    X, y,
    K_outer, K_inner,
    seed,
    model_params
) {

  run_ridge_nested_cv(
    X = X,
    y = y,
    K_outer = K_outer,
    K_inner = K_inner,
    seed = seed,
    standardize = model_params$standardize
  ) |> data.table::as.data.table()
}


