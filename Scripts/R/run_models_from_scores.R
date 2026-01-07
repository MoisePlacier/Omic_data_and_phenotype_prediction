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
    seed = 123
) {
  library(data.table)
  set.seed(seed)

  res_all <- list()

  for (ph in phenotypes) {

    ## ---------- Témoin aléatoire ----------
    snp_random <- sample(colnames(X_full), topK)
    X_rand <- X_full[, snp_random, drop = FALSE]

    perf_rand <- model_fun(
      X = X_rand,
      y = y_full[[ph]],
      K_outer = K_outer,
      K_inner = K_inner,
      seed = seed + 1,
      model_params = model_params
    )

    perf_rand[, `:=`(
      phenotype = ph,
      method = "Random",
      context = "temoin",
      model = model_name,
      topK = topK,
      n_snps = ncol(X_rand)
    )]

    res_all[[length(res_all) + 1]] <- perf_rand

    ## ---------- Méthodes de scoring ----------
    for (m in methods) {
      for (ctx in contexts) {

        dt_sub <- all_scores_dt[phenotype == ph & method == m & context == ctx]

        # Skip si aucune donnée pour cette combinaison
        if (nrow(dt_sub) < 2) next

        snp_sel <- dt_sub[order(rank)][1:min(topK, .N), SNP]

        # Skip si moins de 2 SNP sélectionnés
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
          n_snps = ncol(X_sub)
        )]

        res_all[[length(res_all) + 1]] <- perf
      }
    }
  }

  return(rbindlist(res_all, fill = TRUE))
}

