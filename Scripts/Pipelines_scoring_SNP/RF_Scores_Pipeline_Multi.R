
#!/usr/bin/env Rscript

##########################################
## Commande bash :
## Rscript RF_Scores_Pipeline_Multi.R configs/config_RF_MULTI.yml
##########################################

############################
## Librairies
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
  stop("Usage: Rscript RF_Scores_Pipeline.R config.yml")
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
## Fonction RF
############################
run_RF_analysis <- function(task, X, y, save_dir, base_name) {

  snp_idx   <- task$snp_idx[[1]]
  test_idx  <- task$test_idx[[1]]
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)

  X_train <- X[train_idx, snp_idx, drop = FALSE]
  X_test  <- X[test_idx,  snp_idx, drop = FALSE]
  y_train <- y[train_idx]
  y_test  <- y[test_idx]

  #mtry_val <- min(500, ncol(X_train))

  rf_fit <- ranger(
    x = X_train,
    y = y_train,
    num.trees =  2000,
    mtry = max(floor(ncol(X_train)/50), 500) ,
    min.node.size = 1,
    importance = "impurity_corrected",
    oob.error = FALSE,
    num.threads = 5,
    seed = 123
  )

  y_pred <- predict(rf_fit, data = X_test)$predictions

  perf_dt <- data.table(
    phenotype   = base_name,
    subset_size = task$subset_size,
    iter        = task$iter,
    fold        = task$fold,
    mtry        = mtry_val,
    R2          = cor(y_test, y_pred)^2,
    RMSE        = sqrt(mean((y_test - y_pred)^2))
  )

  imp_dt <- data.table(
    phenotype   = base_name,
    SNP         = colnames(X_train),
    importance  = as.numeric(importance(rf_fit)),
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

  saveRDS(perf_dt, file.path(save_dir, paste0("perf_RF", tag)))
  saveRDS(imp_dt,  file.path(save_dir, paste0("imp_RF",  tag)))

  invisible(NULL)
}

############################
## Boucle sur les phénotypes
############################
for (BASE_NAME in PHENOS) {

  message(">>> Running RF for phenotype: ", BASE_NAME)

  pheno <- pheno_all[, c("ID", BASE_NAME), with = FALSE]

  merged <- merge(pheno, loc[, .(Population, ID)], by = "ID")
  merged <- merge(merged, geno, by = "ID")

  y <- as.numeric(merged[[BASE_NAME]])
  X <- as.matrix(merged[, 23:ncol(merged)])

  SAVE_DIR <- file.path(SAVE_ROOT, paste0("RF_", BASE_NAME))
  dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

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
    seq_len(N_ITER),
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
  N_CORES <- max(2, detectCores() - 2)
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)

  foreach(
    i = seq_len(nrow(tasks)),
    .packages = c("ranger", "data.table", "caret")
  ) %dopar% {
    run_RF_analysis(tasks[i], X, y, SAVE_DIR, BASE_NAME)
  }

  stopCluster(cl)

  ## CONCATÉNATION FINALE

  perf_files <- list.files(
    SAVE_DIR,
    pattern = "^perf_.*\\.rds$",
    full.names = TRUE
  )

  Imp_files <- list.files(
    SAVE_DIR,
    pattern = "^imp_.*\\.rds$",
    full.names = TRUE
  )

  perf_all <- rbindlist(lapply(perf_files, readRDS))
  vip_all  <- rbindlist(lapply(Imp_files,  readRDS))

  saveRDS(
    perf_all,
    file = file.path(SAVE_DIR, paste0("RF_perf_all_", BASE_NAME, ".rds"))
  )

  saveRDS(
    vip_all,
    file = file.path(SAVE_DIR, paste0("RF_scores_all_", BASE_NAME, ".rds"))
  )

  file.remove(perf_files)
  file.remove(Imp_files)
}
