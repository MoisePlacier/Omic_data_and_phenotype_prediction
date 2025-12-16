#!/usr/bin/env Rscript

library(pls)
library(caret)
library(data.table)
library(parallel)
library(doParallel)
library(yaml)

#---------------------------------
# Lecture du YAML
#---------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) stop("Usage: Rscript VIP_Pipeline.R config.yml")

config_file <- args[1]
config <- yaml::read_yaml(config_file)

Y_FILE        <- config$data$phenotype
X_FILE        <- config$data$genomic
LOC_FILE      <- config$data$localisation
BASE_NAME     <- config$analysis$phenotype_name
N_ITER        <- config$analysis$n_iter
SUBSET_SIZE   <- config$analysis$subset_size
SEED          <- config$analysis$seed
K_OUTER       <- config$cv$outer_folds
SAVE_DIR      <- config$output$save_dir

set.seed(SEED)

#---------------------------------
# Chargement des données
#---------------------------------
load(Y_FILE)        # doit créer Phenotype
load(X_FILE)        # doit créer Genomic
load(LOC_FILE)      # doit créer localisation

pheno <- as.data.table(Phenotype)
geno  <- as.data.table(Genomic)
loc   <- as.data.table(localisation)

geno[, ID := rownames(Genomic)]
pheno[, ID := rownames(Phenotype)]

merged <- merge(pheno, loc[, .(Population, ID)], by = "ID")
merged <- merge(merged, geno, by = "ID")

y <- as.numeric(merged[[BASE_NAME]])
X <- as.matrix(merged[, 23:ncol(merged)])

dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

#---------------------------------
# Fonction VIP
#---------------------------------
calc_vip <- function(pls_model) {
  W <- pls_model$loading.weights
  SSY <- colSums(pls_model$Yscores^2)
  p <- nrow(W)
  vip <- numeric(p)
  for(j in 1:p) {
    vip[j] <- sqrt(p * sum(SSY * (W[j,]^2)) / sum(SSY))
  }
  names(vip) <- rownames(W)
  vip
}

run_pls_analysis <- function(task, X, y, save_dir) {
  snp_idx <- task$snp_idx[[1]]
  test_idx <- task$test_idx[[1]]

  X_sub <- X[, snp_idx, drop = FALSE]
  train_idx <- setdiff(seq_len(nrow(X_sub)), test_idx)

  X_train_df <- as.data.frame(X_sub[train_idx, , drop = FALSE])
  y_train <- y[train_idx]
  X_test_df <- as.data.frame(X_sub[test_idx, , drop = FALSE])
  y_test <- y[test_idx]

  names(X_train_df) <- make.names(names(X_train_df))
  names(X_test_df) <- make.names(names(X_test_df))

  inner_train_index <- createDataPartition(y_train, p = 0.9, list = FALSE)
  inner_X_train <- X_train_df[inner_train_index, , drop = FALSE]
  inner_y_train <- y_train[inner_train_index]
  inner_X_test <- X_train_df[-inner_train_index, , drop = FALSE]
  inner_y_test <- y_train[-inner_train_index]

  df_inner_train <- data.frame(y = inner_y_train, inner_X_train)

  PLS2 <- plsr(y ~ ., data = df_inner_train, ncomp = 2, validation = "none")
  PLS3 <- plsr(y ~ ., data = df_inner_train, ncomp = 3, validation = "none")

  y_pred2 <- as.vector(predict(PLS2, newdata = inner_X_test, ncomp = 2))
  y_pred3 <- as.vector(predict(PLS3, newdata = inner_X_test, ncomp = 3))

  RMSE_inner_test_2 <- sqrt(mean((inner_y_test - y_pred2)^2))
  RMSE_inner_test_3 <- sqrt(mean((inner_y_test - y_pred3)^2))

  ncomp_opt <- ifelse(RMSE_inner_test_2 < RMSE_inner_test_3, 2, 3)

  df_train <- data.frame(y = y_train, X_train_df)
  pls_final <- plsr(y ~ ., data = df_train, ncomp = ncomp_opt, validation = "none")

  y_pred_final <- as.vector(predict(pls_final, newdata = X_test_df, ncomp = ncomp_opt))
  R2_test <- cor(y_test, y_pred_final)^2
  RMSE_test <- sqrt(mean((y_test - y_pred_final)^2))

  vip_res <- calc_vip(pls_final)

  perf_dt <- data.table(
    subset_size = task$subset_size, iter = task$iter, fold = task$fold,
    ncomp_opt = ncomp_opt, R2_test = R2_test, RMSE_test = RMSE_test
  )
  vip_dt <- data.table(
    SNP = names(vip_res), VIP = as.numeric(vip_res),
    subset_size = task$subset_size, iter = task$iter, fold = task$fold
  )

  filename_base <- paste0("_s", task$subset_size, "_iter", task$iter, "_fold", task$fold, ".rds")
  saveRDS(perf_dt, file.path(save_dir, paste0("perf", filename_base)))
  saveRDS(vip_dt, file.path(save_dir, paste0("vip", filename_base)))

  invisible(NULL)
}

#---------------------------------
# Préparation des subsets et folds
#---------------------------------
snp_subsets <- replicate(N_ITER, sample(ncol(X), SUBSET_SIZE), simplify = FALSE)
outer_folds <- lapply(seq_len(N_ITER), function(i) createFolds(y, k = K_OUTER, returnTrain = FALSE))

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

#---------------------------------
# Exécution parallèle
#---------------------------------
N_CORES <- max(1, detectCores() - 2)
cl <- makeCluster(N_CORES)
registerDoParallel(cl)

results <- foreach(i = seq_len(nrow(tasks)),
                   .packages = c("pls","caret","data.table"),
                   .export = c("X","y","calc_vip","run_pls_analysis","tasks","SAVE_DIR")) %dopar% {
                     run_pls_analysis(tasks[i], X, y, SAVE_DIR)
                   }

stopCluster(cl)
