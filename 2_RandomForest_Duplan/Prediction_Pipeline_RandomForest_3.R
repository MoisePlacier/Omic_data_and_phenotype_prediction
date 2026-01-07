#!usr/bin/Rscript
# Prediction Pipeline using Random Forest with Grid Search for Hyperparameter Tuning (version 3)

#Define seed
seed(666)

cat("\nStarting execution : \n > Prediction Pipeline using Random Forest with Grid Search for Hyperparameter Tuning (version 3)\n\n\n")

# ---- clean environment ----
rm(list = ls())
while (!is.null(dev.list())) {
  dev.off()
}

# ---- libraries and directories ----

library(anyLib)
pkg <- c("caret", "ranger", "argparse", "mixOmics", "futile.logger", "foreach", "doParallel")
suppressMessages(anyLib(pkg))

home_dir <- file.path("/home", "INRA", "Aduplan")
rf_dir <- file.path(home_dir, "work", "03_prediction", "02_random_forest")
data_dir <- file.path(home_dir, "01_data")
output_dir <- file.path(rf_dir, "output")

# ---- parser ----

parser <- argparse::ArgumentParser(description = "Script to perform a Random Forest training using a grid search for hyperparameter tuning.")

parser$add_argument(
  "--working_dir",
  help = "Directory where you work",
  type = "character",
  default = rf_dir
)

parser$add_argument(
  "-o", "--output_dir",
  help = "Directory where the output files will be saved",
  type = "character",
  default = output_dir
)

parser$add_argument(
  "-d", "--data_path",
  help = "Path to the omics data file",
  type = "character",
  default = file.path(data_dir, "PoplarDataORL_reduced.RData")
)

parser$add_argument(
  "--n_core",
  help = "Number of cores to use for parallel processing",
  type = "integer",
  default = ifelse(Sys.info()["sysname"] == "Linux", 10, 3)
)

parser$add_argument(
  "--fs_proportion",
  help = "Proportion of features to select for the Random Forest model",
  type = "numeric",
  default = 1
)

parser$add_argument(
  "--n_folds",
  help = "Number of folds for cross-validation",
  type = "integer",
  default = 5
)

parser$add_argument(
  "--train_set_size",
  help = "Proportion of the dataset to use for training",
  type = "numeric",
  default = 0.8
)

parser$add_argument(
  "-v", "--verbose",
  action = "count",
  help = "Increase verbosity level (e.g., -v, -vv, -vvv)",
  default = 0
)

args <- parser$parse_args()

# ---- logging ----

flog.layout(layout.format("[~l] ~t - ~m"))
verbosity_levels <- c(ERROR, WARN, INFO, DEBUG)
selected_level <- verbosity_levels[pmin(args$verbose + 1, length(verbosity_levels))]
flog.threshold(selected_level)

# ---- checking if directories exist ----

if (!dir.exists(args$output_dir)) {
  flog.info("Output directory %s does not exist. Creating it.", args$output_dir)
  dir.create(args$output_dir, recursive = TRUE)
}

# ---- source files ----

flog.info("Sourcing scripts from %s", file.path(rf_dir, "R"))
source(file.path(rf_dir, "R", "feature_ranking_spls.R"))
source(file.path(rf_dir, "R", "random-forest_gridsearch.R"))
flog.info("  > Scripts sourced successfully.")

# ---- load data ----

flog.info("Loading omics data from %s", args$data_path)
data <- get(load(args$data_path))
flog.info("  > Data loaded successfully.")

# ---- ensure unique column names ----

prefix <- list(
  "Genomic" = "G_",
  "Transcriptomic" = "T_",
  "Epigenetic_CG" = "ECG_",
  "Epigenetic_CHG" = "ECHG_",
  "Epigenetic_CHH" = "ECHH_",
  "Phenotype" = ""
)

for (group in names(data)) {
  colnames(data[[group]]) <- paste0(prefix[[group]], colnames(data[[group]]))
}

# ---- combine omic data ----

data$Epigenetic <- cbind(
  data$Epigenetic_CG,
  data$Epigenetic_CHG,
  data$Epigenetic_CHH
)
data$Comb_G_T <- cbind(
  data$Genomic,
  data$Transcriptomic
)
data$Comb_E_T <- cbind(
  data$Epigenetic,
  data$Transcriptomic
)
data$Comb_E_G <- cbind(
  data$Epigenetic,
  data$Genomic
)
data$Comb_E_G_T <- cbind(
  data$Genomic,
  data$Transcriptomic,
  data$Epigenetic
)

# ---- skip groups ----

skip_group <- NULL

flog.info("Available omic data groups:")
for (omic in names(data)) {
  flog.info("  > %s", omic)
}
flog.info("Skipping omic data groups: %s", paste(skip_group, collapse = ", "))

# ---- debug ----

flog.debug("Arguments used:")
for (arg in names(args)) {
  flog.debug("  > %s: %s", arg, args[[arg]])
}

# ---- variables initialization ----

cross_val_results <- list()

# ---- prediction ----

traits <- "Phenotype"
if (!traits %in% names(data)) {
  flog.error("Traits '%s' not found in the data.", traits)
  quit(status = 1)
}

global_start_time <- Sys.time()

for (omic in names(data)) {
  if (omic %in% skip_group) {
    flog.info("Skipping omic data: %s", omic)
    next
  }
  if (omic == traits) next
  cross_val_results <- list()
  flog.info("\n-----------------------------------")
  flog.info("Processing omic data: %s", omic)
  for (trait in colnames(data[[traits]])) {
    flog.info("  > Processing trait: %s", trait)
    start_time <- Sys.time()

    X <- data[[omic]]
    Y <- data[[traits]][[trait]]

    param_grid <- expand.grid(
      mtry = round(sqrt(ncol(X))) + c(-5, 0, 5),
      splitrule = c("variance", "extratrees"),
      min.node.size = 5,
      num.trees = c(500, 1000, 1500),
      max.depth = c(50, 65, 75, 90)
    )

    # Check if X and Y are not NULL
    if (is.null(X) || is.null(Y)) {
      flog.warn("  > Skipping trait %s for omic %s due to missing data.", trait, omic)
      next
    }

    grid_search_results <- grid_search_with_fs(
      X = X,
      Y = Y,
      data_proportion = args$fs_proportion,
      param_grid = param_grid,
      parallelize = TRUE,
      train_set_size = args$train_set_size,
      n_fold = args$n_folds,
      n_cores = args$n_core,
      random_seed = 42
    )

    cross_val_results[[trait]] <- list(
      # results = grid_search_results$results,
      best_model = grid_search_results$best_model,
      best_params = grid_search_results$best_params,
      best_rmse = grid_search_results$best_rmse
    )

    end_time <- Sys.time()
    runtime <- format(end_time - start_time, digits = 2, nsmall = 2)
    flog.info("  > Finished processing trait %s for omic %s in %s", trait, omic, runtime)

    # --- save intermediate results (one file per omic) ---

    output_file <- file.path(args$output_dir, sprintf("rf_results_%s.RData", omic))
    flog.info("Saving results to %s", output_file)
    save(cross_val_results, file = output_file)
  }
}

global_end_time <- Sys.time()
runtime_global <- format(global_end_time - global_start_time, digits = 2, nsmall = 2)
flog.info("\n-----\nFinished processing all omic data in %s", runtime_global)
