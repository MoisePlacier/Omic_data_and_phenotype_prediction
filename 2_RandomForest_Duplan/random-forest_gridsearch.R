#' train_rf_grid_search
#'
#' @description Train a Random Forest model with grid search for hyperparameter tuning
#'
#' @param X_train Training features
#' @param Y_train Training response variable
#' @param param_grid Data frame containing hyperparameter combinations to test
#'
#' @return A list containing the results of the grid search, the best model, best parameters, and best RMSE
train_rf_grid_search <- function(X_train, Y_train, param_grid) {
  if (is.null(param_grid$mtry)) {
    param_grid$mtry <- floor(sqrt(ncol(X_train)))
  }
  if (is.null(param_grid$splitrule)) {
    param_grid$splitrule <- "variance"
  }
  if (is.null(param_grid$min.node.size)) {
    param_grid$min.node.size <- 5
  }
  if (is.null(param_grid$num.trees)) {
    param_grid$num.trees <- 500
  }
  if (is.null(param_grid$max.depth)) {
    param_grid$max.depth <- 30
  }

  # Initialize tracking
  results <- list()
  best_rmse <- Inf
  best_model <- NULL
  best_params <- NULL
  i <- 1
  current_step <- 0
  n_steps <- length(unique(param_grid$num.trees)) * length(unique(param_grid$max.depth))
  flog.info("Starting hyperparameter tuning...")
  for (nt in unique(param_grid$num.trees)) {
    for (md in unique(param_grid$max.depth)) {
      current_step <- current_step + 1
      flog.info(
        "[%s/%s] >> Tuning for num.trees = %d, max.depth = %d",
        current_step, n_steps, nt, md
      )

      tunegrid <- unique(param_grid[
        param_grid$num.trees == nt & param_grid$max.depth == md,
        c("mtry", "splitrule", "min.node.size")
      ])

      ctrl <- caret::trainControl(
        method = "cv",
        number = 5,
        savePredictions = "final",
        allowParallel = TRUE
      )

      start_time <- Sys.time()

      rf_model <- caret::train(
        x = X_train,
        y = Y_train,
        method = "ranger",
        trControl = ctrl,
        tuneGrid = tunegrid,
        num.trees = nt,
        max.depth = md,
        importance = "impurity",
        metric = "RMSE"
      )

      end_time <- Sys.time()
      runtime <- format(end_time - start_time, digits = 2, nsmall = 2)
      flog.info("    > runtime : %s", runtime)

      # Compute R² per fold
      fold_r2 <- by(
        data = data.frame(obs = rf_model$pred$obs, pred = rf_model$pred$pred, Resample = rf_model$pred$Resample),
        INDICES = rf_model$pred$Resample,
        FUN = function(df) cor(df$obs, df$pred)^2
      )

      rf_model$fold_r2 <- fold_r2
      rf_model$mean_r2 <- mean(fold_r2)

      results[[i]] <- rf_model
      i <- i + 1

      # Check best model (by RMSE)
      if (min(rf_model$results$RMSE) < best_rmse) {
        best_rmse <- min(rf_model$results$RMSE)
        best_model <- rf_model
        best_params <- list(num.trees = nt, max.depth = md)
      }
    }
  }
  return(list(results = results, best_model = best_model, best_params = best_params, best_rmse = best_rmse))
}

#' train_rf_grid_search_parallel
#'
#' @description Train a Random Forest model with grid search for hyperparameter tuning in parallel
#'
#' @param X_train Training features
#' @param Y_train Training response variable
#' @param param_grid Data frame containing hyperparameter combinations to test
#' @param n_fold Number of folds for cross-validation
#' @param n_cores Number of cores to use for parallel processing (default: all available cores - 1)
#'
#' @return A list containing the results of the grid search, the best model, best parameters, and best RMSE
train_rf_grid_search_parallel <- function(X_train, Y_train,
                                          param_grid = NULL,
                                          n_fold = 5,
                                          n_cores = max(1, parallel::detectCores() - 1)) {
  library(foreach)
  library(doParallel)

  if (n_fold < 2) {
    flog.error("n_fold must be at least 2 for cross-validation.")
    return(NULL)
  }

  # Parameter validation
  if (is.null(param_grid$mtry)) {
    param_grid$mtry <- floor(sqrt(ncol(X_train)))
  }
  if (is.null(param_grid$splitrule)) {
    param_grid$splitrule <- "variance"
  }
  if (is.null(param_grid$min.node.size)) {
    param_grid$min.node.size <- 5
  }
  if (is.null(param_grid$num.trees)) {
    param_grid$num.trees <- 500
  }
  if (is.null(param_grid$max.depth)) {
    param_grid$max.depth <- 30
  }

  # Create all combinations of num.trees and max.depth
  param_combinations <- expand.grid(
    nt = unique(param_grid$num.trees),
    md = unique(param_grid$max.depth)
  )
  n_combinations <- nrow(param_combinations)

  # Set up parallel backend
  if (n_cores == -1 || n_cores >= parallel::detectCores()) {
    flog.info("Using all available cores: %d", parallel::detectCores())
    n_cores <- parallel::detectCores()
  } else if (n_cores <= 0 && n_cores != -1) {
    flog.error(
      "Invalid number of cores specified: %d. Please provide a positive integer or -1 for all available cores.",
      n_cores
    )
    return(NULL)
  }

  # Log system information for debugging
  total_cores <- parallel::detectCores()
  physical_cores <- parallel::detectCores(logical = FALSE)
  flog.info("System info: %d logical cores, %d physical cores", total_cores, physical_cores)
  flog.info("Using %d cores for parallel processing", n_cores)

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Export necessary variables to the cluster
  clusterExport(cl,
    varlist = c("X_train", "Y_train", "param_grid"),
    envir = environment()
  )

  # Load required packages on each node
  clusterEvalQ(cl, {
    library(caret)
    library(ranger)
  })

  flog.info("Starting parallel hyperparameter tuning...")
  global_start_time <- Sys.time()

  # Parallel foreach loop
  results <- foreach(
    i = 1:n_combinations, .packages = c("caret", "ranger"),
    .errorhandling = "pass"
  ) %dopar% {
    nt <- param_combinations$nt[i]
    md <- param_combinations$md[i]

    tunegrid <- unique(param_grid[
      param_grid$num.trees == nt & param_grid$max.depth == md,
      c("mtry", "splitrule", "min.node.size")
    ])

    ctrl <- caret::trainControl(
      method = "cv",
      number = n_fold,
      savePredictions = "final",
      allowParallel = FALSE
    )

    start_time <- Sys.time()

    rf_model <- tryCatch(
      {
        caret::train(
          x = X_train,
          y = Y_train,
          method = "ranger",
          trControl = ctrl,
          tuneGrid = tunegrid,
          num.trees = nt,
          max.depth = md,
          importance = "impurity",
          metric = "RMSE"
        )
      },
      error = function(e) {
        return(list(error = TRUE, message = e$message))
      }
    )

    end_time <- Sys.time()
    runtime <- format(end_time - start_time, digits = 2, nsmall = 2)

    # Check if there was an error
    if (is.list(rf_model) && !is.null(rf_model$error)) {
      return(list(
        success = FALSE,
        params = list(num.trees = nt, max.depth = md),
        error_message = rf_model$message,
        runtime = runtime
      ))
    }

    # Compute R² per fold
    fold_r2 <- by(
      data = data.frame(
        obs = rf_model$pred$obs, pred = rf_model$pred$pred,
        Resample = rf_model$pred$Resample
      ),
      INDICES = rf_model$pred$Resample,
      FUN = function(df) cor(df$obs, df$pred)^2
    )

    rf_model$fold_r2 <- fold_r2
    rf_model$mean_r2 <- mean(fold_r2)

    return(list(
      success = TRUE,
      model = rf_model,
      params = list(num.trees = nt, max.depth = md),
      runtime = runtime
    ))
  }

  global_end_time <- Sys.time()
  global_runtime <- format(global_end_time - global_start_time, digits = 2, nsmall = 2)

  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()

  # Process results
  valid_results <- list()
  best_rmse <- Inf
  best_model <- NULL
  best_params <- NULL

  for (i in 1:length(results)) {
    result <- results[[i]]

    if (result$success) {
      flog.info(
        "[%d/%d] num.trees = %d, max.depth = %d completed in %s seconds",
        i, length(results), result$params$num.trees, result$params$max.depth,
        result$runtime
      )

      model <- result$model
      valid_results[[length(valid_results) + 1]] <- model

      # Check if this is the best model so far
      current_rmse <- min(model$results$RMSE)
      if (current_rmse < best_rmse) {
        best_rmse <- current_rmse
        best_model <- model
        best_params <- result$params
      }
    } else {
      flog.error(
        "[%d/%d] num.trees = %d, max.depth = %d failed: %s",
        i, length(results), result$params$num.trees, result$params$max.depth,
        result$error_message
      )
    }
  }

  flog.info("Grid search completed in %s", global_runtime)

  return(list(
    results = valid_results,
    best_model = best_model,
    best_params = best_params,
    best_rmse = best_rmse
  ))
}

#' grid_search_with_fs
#'
#' @description Perform grid search with feature selection for Random Forest model
#'
#' @param X A matrix or data frame of predictor variables (features).
#' @param Y A vector or matrix of response variables.
#' @param data_proportion Proportion of features to select based on importance (default: 1, no selection).
#' @param param_grid Data frame containing hyperparameter combinations to test.
#' @param parallelize Logical indicating whether to run grid search in parallel (default: TRUE).
#' @param train_set_size Proportion of data to use for training (default: 0.8).
#' @param n_fold Number of folds for cross-validation (default: 5).
#' @param n_cores Number of cores to use for parallel processing (default: all available cores - 1).
#' @param random_seed Optional random seed for reproducibility (default: NULL).
#'
#' @return A list containing the results of the grid search, the best model, best parameters, and best RMSE.
grid_search_with_fs <- function(X, Y, data_proportion = 1,
                                param_grid = NULL,
                                parallelize = TRUE,
                                train_set_size = 0.8,
                                n_fold = 5,
                                n_cores = max(1, parallel::detectCores() - 1),
                                random_seed = NULL) {
  if (is.null(X) || is.null(Y)) {
    flog.error("X or Y is NULL")
    return(NULL)
  }
  if (nrow(X) != length(Y)) {
    flog.error("Number of rows in X (%d) doesn't match length of Y (%d)", nrow(X), length(Y))
    return(NULL)
  }

  if (data_proportion < 1) {
    n_features <- round(data_proportion * ncol(X))
    flog.debug("Selecting %d features based on data_proportion: %f", n_features, data_proportion)
    selected_features <- select_features(X, Y, n_features = n_features)
    if (!all(selected_features %in% colnames(X))) {
      flog.error("Some selected features are not present in the dataset.")
      return(NULL)
    }
    X <- X[, selected_features, drop = FALSE]
  }
  X <- as.data.frame(X)

  if (!is.null(random_seed)) {
    set.seed(42)
  }

  train_index <- createDataPartition(Y, p = train_set_size, list = FALSE)
  X_train <- X[train_index, , drop = FALSE]
  Y_train <- Y[train_index]
  X_test <- X[-train_index, , drop = FALSE]
  Y_test <- Y[-train_index]

  if (parallelize) {
    grid_search_results <- train_rf_grid_search_parallel(X_train, Y_train,
      param_grid = param_grid,
      n_fold = n_fold,
      n_cores = n_cores
    )
  } else {
    grid_search_results <- train_rf_grid_search(X_train, Y_train,
      param_grid = param_grid,
      n_fold = n_fold
    )
  }

  return(grid_search_results)
}
