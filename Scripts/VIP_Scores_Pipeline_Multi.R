#!/usr/bin/env Rscript

##########################################
## Commande bash :
## Rscript VIP_Scores_Pipeline_Multi.R configs/config_VIP_MULTI.yml
##########################################

############################
## Librairies
############################
suppressPackageStartupMessages({
  library(pls)
  library(caret)
  library(data.table)
  library(parallel)
  library(doParallel)
  library(yaml)
})

############################
## Lecture des arguments
############################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript VIP_Scores_Pipeline.R config.yml")
}

config <- yaml::read_yaml(args[1])

############################
## Paramètres YAML
############################
Y_FILE   <- config$data$phenotype
X_FILE   <- config$data$genomic
LOC_FILE <- config$data$localisation
RES_FILE <- config$data$residuals

PHENOS   <- config$analysis$phenotypes
USE_RES  <- config$analysis$use_residuals
N_ITER   <- config$analysis$n_iter
SUBSET_SIZE <- config$analysis$subset_size
SEED     <- config$analysis$seed

K_OUTER  <- config$cv$outer_folds
SAVE_ROOT <- config$output$save_dir

set.seed(SEED)

############################
## Chargement des données
############################
load(Y_FILE)      # Phenotype
load(X_FILE)      # Genomic
load(LOC_FILE)    # localisation

geno <- as.data.table(Genomic)
loc  <- as.data.table(localisation)
geno[, ID := rownames(Genomic)]

if (USE_RES) {
  pheno_all <- as.data.table(readRDS(RES_FILE))
} else {
  pheno_all <- as.data.table(Phenotype)
  pheno_all[, ID := rownames(Phenotype)]
}

############################
## Fonction VIP
############################
calc_vip <- function(pls_model) {
  W   <- pls_model$loading.weights
  SSY <- colSums(pls_model$Yscores^2)
  p   <- nrow(W)

  vip <- numeric(p)
  for (j in seq_len(p)) {
    vip[j] <- sqrt(p * sum(SSY * (W[j, ]^2)) / sum(SSY))
  }
  names(vip) <- rownames(W)
  vip
}

############################
## Fonction PLS + VIP
############################
run_pls_analysis <- function(task, X, y, save_dir, base_name) {

  snp_idx  <- task$snp_idx[[1]]
  test_idx <- task$test_idx[[1]]

  X_sub <- X[, snp_idx, drop = FALSE]
  train_idx <- setdiff(seq_len(nrow(X_sub)), test_idx)

  X_train_df <- as.data.frame(X_sub[train_idx, , drop = FALSE])
  y_train <- y[train_idx]
  X_test_df  <- as.data.frame(X_sub[test_idx, , drop = FALSE])
  y_test <- y[test_idx]

  names(X_train_df) <- make.names(names(X_train_df))
  names(X_test_df)  <- make.names(names(X_test_df))

  ## Sélection simple ncomp (2 vs 3)
  inner_idx <- createDataPartition(y_train, p = 0.9, list = FALSE)
  df_inner  <- data.frame(y = y_train[inner_idx],
                          X_train_df[inner_idx, , drop = FALSE])

  PLS2 <- plsr(y ~ ., data = df_inner, ncomp = 2, validation = "none")
  PLS3 <- plsr(y ~ ., data = df_inner, ncomp = 3, validation = "none")

  y2 <- as.vector(predict(PLS2, X_train_df[-inner_idx, ], ncomp = 2))
  y3 <- as.vector(predict(PLS3, X_train_df[-inner_idx, ], ncomp = 3))

  rmse2 <- sqrt(mean((y_train[-inner_idx] - y2)^2))
  rmse3 <- sqrt(mean((y_train[-inner_idx] - y3)^2))

  ncomp_opt <- ifelse(rmse2 < rmse3, 2, 3)

  ## Modèle final
  df_train <- data.frame(y = y_train, X_train_df)
  pls_final <- plsr(y ~ ., data = df_train, ncomp = ncomp_opt, validation = "none")

  y_pred <- as.vector(predict(pls_final, X_test_df, ncomp = ncomp_opt))

  R2_test   <- cor(y_test, y_pred)^2
  RMSE_test <- sqrt(mean((y_test - y_pred)^2))

  vip_res <- calc_vip(pls_final)

  perf_dt <- data.table(
    phenotype   = base_name,
    subset_size = task$subset_size,
    iter        = task$iter,
    fold        = task$fold,
    ncomp       = ncomp_opt,
    R2          = R2_test,
    RMSE        = RMSE_test
  )

  vip_dt <- data.table(
    phenotype   = base_name,
    SNP         = names(vip_res),
    VIP         = as.numeric(vip_res),
    subset_size = task$subset_size,
    iter        = task$iter,
    fold        = task$fold
  )

  tag <- paste0(
    "_", base_name,
    "_s", task$subset_size,
    "_iter", task$iter,
    "_fold", task$fold,
    ".rds"
  )

  saveRDS(perf_dt, file.path(save_dir, paste0("perf_VIP", tag)))
  saveRDS(vip_dt,  file.path(save_dir, paste0("vip", tag)))

  invisible(NULL)
}

############################
## Boucle phénotypes
############################
for (BASE_NAME in PHENOS) {

  message(">>> VIP – phenotype: ", BASE_NAME)

  pheno <- pheno_all[, c("ID", BASE_NAME), with = FALSE]

  merged <- merge(pheno, loc[, .(Population, ID)], by = "ID")
  merged <- merge(merged, geno, by = "ID")

  y <- as.numeric(merged[[BASE_NAME]])
  X <- as.matrix(merged[, 23:ncol(merged)])

  SAVE_DIR <- file.path(SAVE_ROOT, paste0("VIP_", BASE_NAME))
  dir.create(SAVE_DIR, recursive = TRUE, showWarnings = FALSE)

  snp_subsets <- replicate(
    N_ITER,
    sample(ncol(X), SUBSET_SIZE),
    simplify = FALSE
  )

  outer_folds <- lapply(
    seq_len(N_ITER),
    function(i) createFolds(y, k = K_OUTER, returnTrain = FALSE)
  )

  tasks <- rbindlist(lapply(seq_len(N_ITER), function(i) {
    rbindlist(lapply(seq_len(K_OUTER), function(f) {
      data.table(
        subset_size = SUBSET_SIZE,
        iter = i,
        fold = f,
        snp_idx = list(snp_subsets[[i]]),
        test_idx = list(outer_folds[[i]][[f]])
      )
    }))
  }))

  N_CORES <- max(1, detectCores() - 2)
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)

  foreach(
    i = seq_len(nrow(tasks)),
    .packages = c("pls", "caret", "data.table")
  ) %dopar% {
    run_pls_analysis(tasks[i], X, y, SAVE_DIR, BASE_NAME)
  }

  stopCluster(cl)
}
