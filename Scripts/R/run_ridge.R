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
    standardize = TRUE,
    n_random_repeats = 10  # Nouveau paramètre
) {
  library(data.table)
  cat(">>> ENTERING run_model_from_scores <<<\n")

  set.seed(seed)
  res_all <- list()

  for (ph in phenotypes) {

    ## --- Témoin Aléatoire répété ---
    message("Running Random controls for phenotype: ", ph)
    for (i in seq_len(n_random_repeats)) {
      # On change le seed à chaque répétition pour avoir des tirages différents
      snp_random <- sample(colnames(X_full), topK)
      X_rand <- X_full[, snp_random, drop = FALSE]

      perf_rand <- run_ridge_nested_cv(
        X = X_rand,
        y = y_full,
        K_outer = K_outer,
        K_inner = K_inner,
        seed = seed + i,
        standardize = standardize
      )

      perf_rand_dt <- as.data.table(perf_rand)
      perf_rand_dt[, `:=`(
        phenotype = ph,
        method = "Random",
        context = paste0("repeat_", i), # Pour identifier les différents tirages
        topK = topK,
        n_snps = ncol(X_rand)
      )]
      res_all[[length(res_all) + 1]] <- perf_rand_dt
    }

    ## --- Boucle sur tes méthodes expertes ---
    for (m in methods) {
      for (ctx in contexts) {

        message("Running ridge | ph = ", ph, " | method = ", m, " | ctx = ", ctx)

        # Sélection des TOP K SNPs basée sur le rang
        snp_sel <- all_scores_dt[
          phenotype == ph & method == m & context == ctx
        ][order(rank)][1:topK, SNP]

        # Sécurité : vérifier que les SNPs sont bien dans la matrice de génotypage
        snp_sel <- snp_sel[snp_sel %in% colnames(X_full)]

        if (length(snp_sel) < 2) {
          warning("Pas assez de SNPs trouvés pour ", m, " ", ctx)
          next
        }

        X_sub <- X_full[, snp_sel, drop = FALSE]

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




run_hybrid_pipeline <- function(
    all_scores_dt,
    phenotypes,
    X_full,
    y_full,
    topK = 500,
    method_gwas = "FarmCPU", context_gwas = "GWAS",
    method_rf = "RF", context_rf = "LMM_residuals",
    K_outer = 5, K_inner = 5,
    seed = 123
) {
  library(data.table)
  set.seed(seed)
  res_hybrid <- list()

  for (ph in phenotypes) {
    message(">>> Processing Hybrid Model for: ", ph)

    # 1. Sélection dirigée : GWAS + RF
    pool_gwas <- all_scores_dt[phenotype == ph & method == method_gwas & context == context_gwas][order(rank)]
    pool_rf   <- all_scores_dt[phenotype == ph & method == method_rf & context == context_rf][order(rank)]

    if(nrow(pool_gwas) == 0 | nrow(pool_rf) == 0) {
      warning("Scores manquants pour le phénotype: ", ph)
      next
    }

    # On prend la moitié de chaque
    half <- floor(topK / 2)
    snps_gwas <- pool_gwas[1:half, SNP]
    snps_rf   <- pool_rf[1:half, SNP]

    # Fusion et complétion pour atteindre exactement topK unique
    combined_snps <- unique(c(snps_gwas, snps_rf))

    if (length(combined_snps) < topK) {
      needed <- topK - length(combined_snps)
      # On complète avec les suivants du GWAS (ou RF) pour boucher les trous
      extra <- setdiff(pool_gwas[(half+1):(half+needed+50), SNP], combined_snps)[1:needed]
      combined_snps <- c(combined_snps, extra)
    }

    # Filtrage sécurité matrice
    combined_snps <- combined_snps[combined_snps %in% colnames(X_full)]

    # 2. Run du modèle Hybride
    perf_hyb <- run_ridge_nested_cv(
      X = X_full[, combined_snps, drop = FALSE],
      y = y_full[[ph]],
      K_outer = K_outer, K_inner = K_inner,
      seed = seed
    )

    perf_hyb_dt <- as.data.table(perf_hyb)
    perf_hyb_dt[, `:=`(
      phenotype = ph,
      method = "Hybrid_Fusion",
      context = paste0(method_gwas, "+", method_rf),
      model = "ridge",
      topK = topK,
      n_snps = length(combined_snps)
    )]
    res_hybrid[[length(res_hybrid) + 1]] <- perf_hyb_dt

    # 3. Témoin aléatoire dédié (même taille que l'hybride)
    snp_rand <- sample(colnames(X_full), length(combined_snps))
    perf_rand <- run_ridge_nested_cv(
      X = X_full[, snp_rand, drop = FALSE],
      y = y_full[[ph]],
      K_outer = K_outer, K_inner = K_inner,
      seed = seed + 1
    )

    perf_rand_dt <- as.data.table(perf_rand)
    perf_rand_dt[, `:=`(
      phenotype = ph,
      method = "Random_Hybrid",
      context = "temoin",
      model = "ridge",
      topK = topK,
      n_snps = length(snp_rand)
    )]
    res_hybrid[[length(res_hybrid) + 1]] <- perf_rand_dt
  }

  return(rbindlist(res_hybrid))
}







################################################################################
# FONCTION : run_elbow_pipeline (Version Répétée)
# Itère sur K et répète la Cross-Validation pour stabiliser les boxplots
################################################################################
run_elbow_pipeline_ridge <- function(
    topK_vector = c(10, 20, 50, 100, 200, 500, 1000, 2000),
    all_scores_dt,
    phenotypes,
    methods,
    contexts,
    X_full,
    y_full,
    K_outer = 5,
    K_inner = 5,
    n_cv_repeats = 10,  # Nombre de fois qu'on refait la 5-fold CV
    seed = 123,
    n_random_repeats = 5
) {
  library(data.table)
  all_elbow_results <- list()

  for (k_val in topK_vector) {
    message("\n--- Calcul pour K = ", k_val, " (", n_cv_repeats, " répétitions de CV) ---")

    # Liste temporaire pour stocker les répétitions de ce palier K
    list_reps_k <- list()

    for (r in seq_len(n_cv_repeats)) {
      # On décale le seed à chaque répétition pour que 'createFolds'
      # à l'intérieur de run_ridge_nested_cv génère des splits différents.
      current_seed <- seed + (r * 1000)

      res_k_rep <- run_ridge_from_scores(
        all_scores_dt = all_scores_dt,
        phenotypes = phenotypes,
        methods = methods,
        contexts = contexts,
        topK = k_val,
        X_full = X_full,
        y_full = y_full,
        K_outer = K_outer,
        K_inner = K_inner,
        seed = current_seed,
        n_random_repeats = n_random_repeats
      )

      # On marque le numéro de la répétition
      res_k_rep[, cv_repeat_id := r]
      list_reps_k[[r]] <- res_k_rep
    }

    # Fusion des répétitions pour ce palier K
    all_elbow_results[[as.character(k_val)]] <- rbindlist(list_reps_k)
  }

  # Fusion finale de tous les paliers K
  final_dt <- rbindlist(all_elbow_results)

  # On s'assure que topK est un facteur pour un bel affichage en boxplot
  final_dt[, topK_factor := factor(topK, levels = sort(unique(topK)))]

  # Nettoyage des témoins Random (comme discuté, on les concatène tous)
  final_dt[method %like% "Random", `:=`(method = "Random", context = "Baseline")]

  return(final_dt)
}
