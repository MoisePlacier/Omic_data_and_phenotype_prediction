library(caret)
library(ranger)
library(data.table)
library(parallel)
library(doParallel)
library(glmnet)

run_model_from_scores <- function(
    all_scores_dt,
    phenotypes,
    methods,
    contexts,
    topK,
    X_full,
    y_full,
    model_name,
    model_fun,
    model_params,
    K_outer = 5,
    K_inner = 5,
    seed = 123,
    n_random_repeats = 5
) {
  library(data.table)
  set.seed(seed)
  res_all <- list()

  for (ph in phenotypes) {
    ## ---------- Témoins aléatoires (Répétés) ----------
    for (i in seq_len(n_random_repeats)) {
      snp_random <- sample(colnames(X_full), topK)
      X_rand <- X_full[, snp_random, drop = FALSE]

      perf_rand <- model_fun(
        X = X_rand,
        y = y_full[[ph]],
        K_outer = K_outer,
        K_inner = K_inner,
        seed = seed + i, # Seed différente pour chaque tirage aléatoire
        model_params = model_params
      )

      perf_rand[, `:=`(
        phenotype = ph,
        method = "Random",
        context = "temoin",
        model = model_name,
        topK = topK,
        n_snps = ncol(X_rand),
        random_rep = i
      )]
      res_all[[length(res_all) + 1]] <- perf_rand
    }

    ## ---------- Méthodes de scoring ----------
    for (m in methods) {
      for (ctx in contexts) {
        dt_sub <- all_scores_dt[phenotype == ph & method == m & context == ctx]
        if (nrow(dt_sub) < 2) next

        snp_sel <- dt_sub[order(rank)][1:min(topK, .N), SNP]
        if (length(snp_sel) < 2) next

        X_sub <- X_full[, colnames(X_full) %in% snp_sel, drop = FALSE]

        perf <- model_fun(
          X = X_sub,
          y = y_full[[ph]],
          K_outer = K_outer,
          K_inner = K_inner,
          seed = seed,
          model_params = model_params
        )

        perf[, `:=`(
          phenotype = ph,
          method = m,
          context = ctx,
          model = model_name,
          topK = topK,
          n_snps = ncol(X_sub),
          random_rep = NA
        )]
        res_all[[length(res_all) + 1]] <- perf
      }
    }
  }
  return(rbindlist(res_all, fill = TRUE))
}



run_elbow_pipeline <- function(
    topK_vector,
    all_scores_dt,
    phenotypes,
    methods,
    contexts,
    X_full,
    y_full,
    model_name = "RF",
    model_fun = run_rf_model,
    model_params = list(rf_grid = rf_grid),
    n_cv_repeats = 5,
    n_random_repeats = 3,
    K_outer = 5,
    K_inner = 5,
    seed = 123
) {
  all_elbow_results <- list()

  for (k_val in topK_vector) {
    message("\n>>> Calcul Palier K = ", k_val)
    list_reps <- list()

    for (r in seq_len(n_cv_repeats)) {
      res_k_rep <- run_model_from_scores(
        all_scores_dt = all_scores_dt,
        phenotypes = phenotypes,
        methods = methods,
        contexts = contexts,
        topK = k_val,
        X_full = X_full,
        y_full = y_full,
        model_name = model_name,
        model_fun = model_fun,
        model_params = model_params,
        K_outer = K_outer,
        K_inner = K_inner,
        seed = seed + (r * 1000),
        n_random_repeats = n_random_repeats
      )
      res_k_rep[, cv_repeat_id := r]
      list_reps[[r]] <- res_k_rep
    }
    all_elbow_results[[as.character(k_val)]] <- rbindlist(list_reps)
  }

  final_dt <- rbindlist(all_elbow_results, fill = TRUE)
  return(final_dt)
}
