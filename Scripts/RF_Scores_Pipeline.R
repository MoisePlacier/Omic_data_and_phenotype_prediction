#!/usr/bin/env Rscript

##########################################
## Commande bash pour run le pipeline :
# Rscript RF_Scores_Pipeline.R configs/config_RF_CIRC2009.yml
##########################################
############################
## Chargement des librairies
############################
suppressPackageStartupMessages({
  library(data.table)
  library(caret)
  library(ranger)
  library(parallel)
  library(doParallel)
  library(yaml)
})

############################
## Lecture des arguments
############################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript run_RF_pipeline.R config_RF.yml")
}

config <- yaml::read_yaml(args[1])

Y_FILE      <- config$data$phenotype
X_FILE      <- config$data$genomic
LOC_FILE    <- config$data$localisation

BASE_NAME   <- config$analysis$phenotype_name
N_ITER      <- config$analysis$n_iter
SUBSET_SIZE <- config$analysis$subset_size
SEED        <- config$analysis$seed

K_OUTER     <- config$cv$outer_folds

SAVE_DIR    <- config$output$save_dir

set.seed(SEED)


############################
## Chargement des données
############################
load(Y_FILE)        # doit créer Phenotype
load(X_FILE)        # doit créer Genomic
load(LOC_FILE)      # doit créer localisation

readRDS("results/Residuals/LM_Pheno_pop_Residuals.rds")

#################### pour faire sur les résiduals :
#pheno <- as.data.table(readRDS("results/Residuals/LM_Pheno_pop_Residuals.rds"))
#################### pour faire sur phénotype !
pheno  <- as.data.table(Phenotype)
pheno[, ID := rownames(Phenotype)]
######################################

geno  <- as.data.table(Genomic)
loc   <- as.data.table(localisation)

geno[, ID := rownames(Genomic)]



merged <- merge(pheno, loc[, .(Population, ID)], by = "ID")
merged <- merge(merged, geno, by = "ID")

y <- as.numeric(merged[[BASE_NAME]])
X <- as.matrix(merged[, 23:ncol(merged)])

############################
## Y et X
############################
y <- as.numeric(merged[[BASE_NAME]])
X <- as.matrix(merged[, 23:ncol(merged)])

############################
## Création dossier sortie
############################

dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

############################
## Fonction RF
############################
run_RF_analysis <- function(task, X, y, save_dir, base_name) {

  snp_idx  <- task$snp_idx[[1]]
  test_idx <- task$test_idx[[1]]
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)

  X_train <- X[train_idx, snp_idx, drop = FALSE]
  X_test  <- X[test_idx,  snp_idx, drop = FALSE]
  y_train <- y[train_idx]
  y_test  <- y[test_idx]

  mtry_val <- min(500, ncol(X_train))

  rf_fit <- ranger(
    x = X_train,
    y = y_train,
    num.trees = 500,
    mtry = mtry_val,
    min.node.size = 5,
    importance = "impurity_corrected",
    oob.error = FALSE,
    num.threads = 4,
    seed = 123
  )

  y_pred <- predict(rf_fit, data = X_test)$predictions

  perf_dt <- data.table(
    phenotype = base_name,
    subset_size = task$subset_size,
    iter = task$iter,
    fold = task$fold,
    mtry = mtry_val,
    R2 = cor(y_test, y_pred)^2,
    RMSE = sqrt(mean((y_test - y_pred)^2))
  )

  imp_dt <- data.table(
    phenotype = base_name,
    SNP = colnames(X_train),
    importance = as.numeric(importance(rf_fit)),
    subset_size = task$subset_size,
    iter = task$iter,
    fold = task$fold
  )

  tag <- paste0(
    "_", base_name,
    "_s", task$subset_size,
    "_iter", task$iter,
    "_fold", task$fold,
    ".rds"
  )

  saveRDS(perf_dt, file.path(save_dir, paste0("perf_RF", tag)))
  saveRDS(imp_dt,  file.path(save_dir, paste0("imp_RF",  tag)))

  invisible(NULL)
}


N_CORES  <- max(1, detectCores() - 2)

############################
## Génération des subsets SNP
############################
snp_subsets <- replicate(
  N_ITER,
  sample(ncol(X), SUBSET_SIZE),
  simplify = FALSE
)

############################
## Folds externes
############################
outer_folds <- lapply(
  1:N_ITER,
  function(i) createFolds(y, k = K_OUTER, returnTrain = FALSE)
)

############################
## Table des tâches
############################
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

############################
## Parallélisation
############################
cl <- makeCluster(N_CORES)
registerDoParallel(cl)

foreach(
  i = seq_len(nrow(tasks)),
  .packages = c("ranger", "data.table", "caret")
) %dopar% {
  run_RF_analysis(tasks[i], X, y, SAVE_DIR, BASE_NAME)
}

stopCluster(cl)
