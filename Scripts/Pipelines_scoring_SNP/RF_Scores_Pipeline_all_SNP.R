############################
## Lecture des arguments
############################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript RF_Scores_Pipeline.R config.yml")
}

library(data.table)
library(ranger)
library(yaml)

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
SEED     <- config$analysis$seed

SAVE_ROOT <- config$output$save_dir

set.seed(SEED)

############################
## Chargement des données
############################
load(Y_FILE)      # Phenotype
load(X_FILE)      # Genomic
load(LOC_FILE)    # localisation

geno <- as.data.table(Genomic)
geno[, ID := rownames(Genomic)]

loc <- as.data.table(localisation)

if (USE_RES) {
  pheno_all <- as.data.table(readRDS(RES_FILE))
} else {
  pheno_all <- as.data.table(Phenotype)
  pheno_all[, ID := rownames(Phenotype)]
}

############################
## Fonction RF importance
############################
run_RF_importance <- function(X, y, phenotype, save_dir) {

  p <- ncol(X)

  rf_fit <- ranger(
    x = X,
    y = y,
    num.trees = 3000,
    mtry = 1000,
    min.node.size = 50,
    max.depth = 6,
    sample.fraction = 0.7,
    importance = "impurity",
    replace = TRUE,
    oob.error = FALSE
  )

  imp_dt <- data.table(
    phenotype  = phenotype,
    SNP        = colnames(X),
    importance = as.numeric(importance(rf_fit))
  )

  saveRDS(
    imp_dt,
    file = file.path(save_dir, paste0("RF_importance_", phenotype, ".rds"))
  )

  invisible(NULL)
}

############################
## Boucle sur les phénotypes
############################
for (BASE_NAME in PHENOS) {

  message(">>> RF importance for phenotype: ", BASE_NAME)

  pheno <- pheno_all[, .(ID, value = get(BASE_NAME))]

  merged <- merge(pheno, loc[, .(Population, ID)], by = "ID")
  merged <- merge(merged, geno, by = "ID")

  y <- as.numeric(merged$value)
  X <- as.matrix(merged[, 23:ncol(merged)])

  SAVE_DIR <- file.path(SAVE_ROOT, paste0("RF_", BASE_NAME))
  dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

  run_RF_importance(X, y, BASE_NAME, SAVE_DIR)
}

