
library(caret)
library(ranger)
library(data.table)


run_rf_nested_cv <- function(
    X, y,
    rf_grid,
    K_outer = 5, K_inner = 5,
    seed = 123
) {
  library(data.table)
  library(caret)
  library(ranger)

  set.seed(seed)
  outer_folds <- createFolds(y, k = K_outer, returnTrain = TRUE)

  results <- data.table()
  predictions_all <- data.table()

  for (f in seq_len(K_outer)) {
    train_idx <- outer_folds[[f]]
    test_idx  <- setdiff(seq_along(y), train_idx)

    X_train <- as.data.frame(X[train_idx, , drop = FALSE])
    X_test  <- as.data.frame(X[test_idx, , drop = FALSE])
    y_train <- y[train_idx]
    y_test  <- y[test_idx]

    # --- CV Interne  ---
    inner_folds <- createFolds(y_train, k = K_inner, returnTrain = TRUE)
    grid_perf <- as.data.table(copy(rf_grid))
    grid_perf[, RMSE_cv := NA_real_]

    for (g in seq_len(nrow(grid_perf))) {
      rmse_inner <- numeric(K_inner)
      for (k in seq_len(K_inner)) {
        tr_idx <- inner_folds[[k]]
        te_idx <- setdiff(seq_along(y_train), tr_idx)

        rf_inner <- ranger(
          x = X_train[tr_idx, ], y = y_train[tr_idx],
          num.trees = grid_perf$num.trees[g],
          mtry = min(grid_perf$mtry[g], ncol(X_train)),
          min.node.size = grid_perf$min.node.size[g],
          num.threads = 3
        )
        y_pred_inner <- predict(rf_inner, X_train[te_idx, ])$predictions
        rmse_inner[k] <- sqrt(mean((y_train[te_idx] - y_pred_inner)^2))
      }
      grid_perf$RMSE_cv[g] <- mean(rmse_inner)
    }

    best <- grid_perf[which.min(RMSE_cv)]

    # --- Modèle Final (Outer Fold) ---
    rf_final <- ranger(
      x = X_train, y = y_train,
      num.trees = best$num.trees,
      mtry = min(best$mtry, ncol(X_train)),
      min.node.size = best$min.node.size,
      importance = "permutation",
      num.threads = 3
    )

    y_pred <- predict(rf_final, X_test)$predictions

    # Calcul des métriques
    R2_test <- 1 - mean((y_test - y_pred)^2) / var(y_test)
    RMSE_test <- sqrt(mean((y_test - y_pred)^2))

    results <- rbind(results, data.table(
      outer_fold = f,
      R2 = R2_test,
      RMSE = RMSE_test,
      num.trees = best$num.trees,
      mtry = best$mtry,
      min.node.size = best$min.node.size
    ))
  }

  return(list(performance = results, predictions = predictions_all))
}



################################################################################

run_rf_from_scores <- function(
    all_scores_dt,
    phenotypes,
    methods,
    contexts,
    topK,
    X_full,
    y_full,
    rf_grid,
    K_outer = 5,
    K_inner = 5,
    seed = 123
) {

  library(data.table)
  set.seed(seed)
  rf_grid <- data.table(
    num.trees = c(500, 1000),
    mtry = c(10, 50),
    min.node.size = c(1, 10)
  )


  res_all <- list()

  for (ph in phenotypes) {
    print(cat("running pheno :", ph, "\n"))

    ## ---------- Témoin aléatoire ----------
    snp_random <- sample(colnames(X_full), topK)
    X_rand <- X_full[, snp_random, drop = FALSE]

    perf_rand <- run_rf_nested_cv(
      X = X_rand,
      y = y_full,
      rf_grid = rf_grid,
      K_outer = K_outer,
      K_inner = K_inner,
      seed = seed + 1
    )

    perf_rand[, `:=`(
      phenotype = ph,
      method = "Random",
      context = "temoin",
      topK = topK,
      n_snps = ncol(X_rand)
    )]

    res_all[[length(res_all) + 1]] <- perf_rand

    for (m in methods) {
      for (ctx in contexts) {

        message(
          "RF | phenotype = ", ph,
          " | method = ", m,
          " | context = ", ctx
        )

        cat(
          "RF | phenotype = ", ph,
          " | method = ", m,
          " | context = ", ctx
        )

        snp_sel <- all_scores_dt[
          phenotype == ph &
            method == m &
            context == ctx
        ][order(rank)][1:topK, SNP]

        if (length(snp_sel) < 2) next

        X_sub <- X_full[, colnames(X_full) %in% snp_sel, drop = FALSE]

        perf <- run_rf_nested_cv(
          X = X_sub,
          y = y_full,
          rf_grid = rf_grid,
          K_outer = K_outer,
          K_inner = K_inner,
          seed = seed
        )

        perf[, `:=`(
          phenotype = ph,
          method = m,
          context = ctx,
          topK = topK,
          n_snps = ncol(X_sub)
        )]

        res_all[[length(res_all) + 1]] <- perf
      }
    }
  }

  return(rbindlist(res_all, fill = TRUE))
}

################################################################################

run_rf_model <- function(
    X, y,
    K_outer, K_inner,
    seed,
    model_params
) {

  out <- run_rf_nested_cv(
    X = X,
    y = y,
    rf_grid = model_params$rf_grid,
    K_outer = K_outer,
    K_inner = K_inner,
    seed = seed
  )

  data.table::as.data.table(out$performance)
}


rf_grid <- data.table(
  num.trees = c(500, 1000),
  mtry = c(10, 50),
  min.node.size = c(1, 10)
)





run_rf_hybrid_pipeline <- function(
    all_scores_dt,
    phenotypes,
    X_full,
    y_full,      # Liste de phénotypes comme dans ton pipeline
    topK = 500,
    rf_grid,
    method_gwas = "FarmCPU", context_gwas = "GWAS",
    method_rf = "RF", context_rf = "LMM_residuals",
    K_outer = 5, K_inner = 5,
    seed = 123,
    n_random_repeats = 5 # RF est lent, on peut réduire un peu les répétitions
) {
  library(data.table)
  set.seed(seed)
  res_hybrid <- list()

  for (ph in phenotypes) {
    message("\n>>> RF Hybrid Pipeline | Phenotype: ", ph)

    # 1. Sélection des deux pools de SNPs
    pool_gwas <- all_scores_dt[phenotype == ph & method == method_gwas & context == context_gwas][order(rank)]
    pool_rf   <- all_scores_dt[phenotype == ph & method == method_rf & context == context_rf][order(rank)]

    if(nrow(pool_gwas) == 0 | nrow(pool_rf) == 0) {
      warning("Scores manquants pour : ", ph)
      next
    }

    # 2. Construction du set Hybride (Top K unique)
    half <- floor(topK / 2)
    snps_gwas <- pool_gwas[1:half, SNP]
    snps_rf   <- pool_rf[1:half, SNP]

    combined_snps <- unique(c(snps_gwas, snps_rf))

    # Complétion si nécessaire pour atteindre exactement topK
    if (length(combined_snps) < topK) {
      needed <- topK - length(combined_snps)
      extra <- setdiff(pool_gwas[(half+1):(half+needed+100), SNP], combined_snps)[1:needed]
      combined_snps <- c(combined_snps, extra)
    }

    combined_snps <- combined_snps[combined_snps %in% colnames(X_full)]

    # 3. Entraînement RF Hybride
    message("... Training Hybrid RF Model")
    perf_hyb <- run_rf_model(
      X = X_full[, combined_snps, drop = FALSE],
      y = y_full[[ph]],
      K_outer = K_outer, K_inner = K_inner,
      seed = seed,
      model_params = list(rf_grid = rf_grid)
    )

    perf_hyb[, `:=`(
      phenotype = ph,
      method = "Hybrid_Fusion",
      context = paste0(method_gwas, "+", method_rf),
      model = "RF",
      topK = topK,
      n_snps = length(combined_snps)
    )]
    res_hybrid[[length(res_hybrid) + 1]] <- perf_hyb

    # 4. Témoins Aléatoires (Répétés)
    message("... Training Random RF Controls (", n_random_repeats, " repeats)")
    for (i in seq_len(n_random_repeats)) {
      snp_rand <- sample(colnames(X_full), length(combined_snps))

      perf_rand <- run_rf_model(
        X = X_full[, snp_rand, drop = FALSE],
        y = y_full[[ph]],
        K_outer = K_outer, K_inner = K_inner,
        seed = seed + i,
        model_params = list(rf_grid = rf_grid)
      )

      perf_rand[, `:=`(
        phenotype = ph,
        method = "Random_Hybrid",
        context = "temoin",
        model = "RF",
        topK = topK,
        n_snps = length(snp_rand)
      )]
      res_hybrid[[length(res_hybrid) + 1]] <- perf_rand
    }
  }

  return(rbindlist(res_hybrid))
}



################################################################################
# FONCTION : run_elbow_pipeline_RF
# Adaptée pour le pipeline Random Forest avec Repeated Nested CV
################################################################################
run_elbow_pipeline_RF <- function(
    topK_vector = c(10, 20, 50, 100, 200, 500), # RF est lent, on limite souvent les K au début
    all_scores_dt,
    phenotypes,
    methods,
    contexts,
    X_full,
    y_full,      # Liste de phénotypes
    rf_grid,     # La grille d'hyperparamètres
    K_outer = 5,
    K_inner = 5,
    n_cv_repeats = 5, # Nombre de répétitions de la 5-fold CV
    seed = 123,
    n_random_repeats = 3 # Nombre de témoins aléatoires par palier K
) {
  library(data.table)
  all_elbow_results <- list()

  for (k_val in topK_vector) {
    message("\n" , paste0(rep("=", 50), collapse = ""))
    message(">>> RF ELBOW PALIER K = ", k_val, " | ", n_cv_repeats, " CV repeats")
    message(paste0(rep("=", 50), collapse = ""))

    list_reps_k <- list()

    for (r in seq_len(n_cv_repeats)) {
      current_seed <- seed + (r * 1000)

      # Appel de ta fonction RF
      # Note : Ta fonction run_rf_from_scores actuelle renvoie un data.table
      # mais attention, à l'intérieur elle doit extraire $performance du nested_cv
      res_k_rep <- run_rf_from_scores(
        all_scores_dt = all_scores_dt,
        phenotypes = phenotypes,
        methods = methods,
        contexts = contexts,
        topK = k_val,
        X_full = X_full,
        y_full = y_full[[phenotypes[1]]], # Adaptation si tu passes un seul phéno
        rf_grid = rf_grid,
        K_outer = K_outer,
        K_inner = K_inner,
        seed = current_seed
      )

      res_k_rep[, cv_repeat_id := r]
      list_reps_k[[r]] <- res_k_rep
    }

    all_elbow_results[[as.character(k_val)]] <- rbindlist(list_reps_k)
  }

  # Fusion finale
  final_dt <- rbindlist(all_elbow_results, fill = TRUE)

  # Création du facteur ordonné pour le plot
  final_dt[, topK := as.numeric(as.character(topK))]
  final_dt[, topK_factor := factor(topK, levels = sort(unique(topK)))]

  # Unification du Random pour le plot
  final_dt[method %like% "Random", `:=`(method = "Random", context = "Baseline")]

  return(final_dt)
}
