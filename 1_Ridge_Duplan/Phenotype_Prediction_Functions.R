# Abdou Rahmane Wade
# abdou.wade@inrae.fr
# INRAE, UMR BioForA, UMR 0588, F-45075 Orleans, France
# 2021

# Phenotype prediction functions  
# For the article :
  # eQTLs are key players in the integration of 
  # genomic and transcriptomic data for phenotype prediction


# Function to print the start of prediction calcul  -------------------------------
#' @descriptionFunction to print the start of prediction calcul
#' @param First.Omic phenotype vector
#' @param Second.Omic number of outer fold of the Nested Cross Validation
#' @param Third.Omic number of inner fold of the Nested Cross Validation
print_message <- function(method_names, trait, First.Omic, Second.Omic = NULL, Third.Omic = NULL) {
  # Construire le message en fonction des valeurs de Second.Omic et Third.Omic
  message <- paste("Start", method_names, "Phenotypic prediction for",trait, "using", First.Omic)
  
  if (Second.Omic != "NULL") {
    message <- paste(message, "and", Second.Omic)
  }
  
  if (Third.Omic != "NULL") {
    message <- paste(message, "and", Third.Omic)
  }
  
  print(message)
}

# Fold Sampling for Nested Cross Validation -------------------------------
#' @description Fold Sampling for Nested Cross Validation
#' @param Y phenotype vector
#' @param outer.fold number of outer fold of the Nested Cross Validation
#' @param inner.fold number of inner fold of the Nested Cross Validation
#' @param iter number of iteration of the Nested Cross Validation
#' @return a list of inner sample for each outer fold of each iteration
Ncv_Fold_Sampling <- function(Y, outer.fold = 5,inner.fold=10, iter = 50){
  sel <- !is.na(Y)
  Y <- Y[sel]
  
  fold <- list()
  for(i in 1:iter)
  {
    
    outer.foldid <- sample(1:outer.fold,size=length(Y),replace=TRUE, 
                           prob = rep(1/outer.fold,outer.fold))
    in.id <- list()
    
    for(outer.t in 1:outer.fold)
    {
      
      outer.train <- which(outer.foldid!=outer.t)
      Y.outer.train <- Y[outer.train]
      in.id[[outer.t]]=sample(1:inner.fold,size=length(Y.outer.train),replace=TRUE, 
                              prob = rep(1/inner.fold,inner.fold))
      
    }
    
    fold[[i]] <- list("out.id"= outer.foldid,
                      "in.id" = in.id)
    
  }
  print("Nested Cross Validation Fold Sampling...Done")
  return(fold)
}
Ncv_RidgeReg <- function(Y, First.Omic, Second.Omic=NULL, cores = 1, alpha = 0, CV.fold){
  # Package Requirements
  require(doSNOW)
  require(glmnet)
  require(caret)
  
  #checks
  X <- First.Omic
  Z <- Second.Omic
  stopifnot(length(Y) == nrow(X))
  if(!is.null(Z)){
    stopifnot(length(Y) == nrow(Z))
  }
  sel <- !is.na(Y)
  X <- X[sel,]
  Y <- Y[sel]
  if(!is.null(Z)){
    Z <- Z[sel,]
  }
  
  # Nested cross-validation parameters
  iter <- length(CV.fold)
  outer.fold <- length(unique(CV.fold[[1]]$out.id))
  
  # Set parallel computing
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  # Computing Regression
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("glmnet","caret"), .options.snow=opts) %dopar% {
    
    # Output object
    beta <- list(); lambda.best <- mu <- R2 <- RMSE <- Rho <- VarA <- VarE <- numeric()
    
    # Outer-fold
    outer.foldid=CV.fold[[i]]$out.id
    
    # Outer-fold loop
    for(outer.t in 1:outer.fold){
      
      outer.test <- which(outer.foldid==outer.t)
      X.outer.test <- X[outer.test,]
      Y.outer.test <- Y[outer.test]
      outer.train <- which(outer.foldid!=outer.t)
      Y.outer.train <- Y[outer.train]
      X.outer.train <- X[outer.train,]
      if(!is.null(Z)){
        Z.outer.test <- Z[outer.test,]
        Z.outer.train <- Z[outer.train,]
      }
      
      # Inner-fold
      inner.foldid=CV.fold[[i]]$in.id[[outer.t]]
      
      # Fit inner model using inner-fold CV: model
      if(is.null(Z)){
        CV.inner.model <- cv.glmnet(x=X.outer.train,y=Y.outer.train,foldid=inner.foldid,alpha=alpha,
                                    keep = TRUE)
        
        lambda.best[outer.t] <- CV.inner.model$lambda.min
        beta[[outer.t]] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[-1]
        mu[outer.t] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[1]
        R2[outer.t] <- R2(mu[outer.t] + (X.outer.test%*%beta[[outer.t]]) , Y.outer.test)
        RMSE[outer.t] <- RMSE(mu[outer.t] + (X.outer.test%*%beta[[outer.t]]) , Y.outer.test)
        Rho[outer.t] <- cor(mu[outer.t] + (X.outer.test%*%beta[[outer.t]]) , Y.outer.test)
        VarA[outer.t] <- var(X.outer.test%*%beta[[outer.t]])
        VarE[outer.t] <- var(Y.outer.train)-VarA[outer.t]
        
      }else{
        
        # Error Msg
        if(is.null(Z)){
          stop("Need Second omic data")
        }
        
        # Fit inner model using inner-fold CV: model
        CV.inner.model <- cv.glmnet(x=cbind(X.outer.train,Z.outer.train),
                                    y=Y.outer.train,foldid=inner.foldid,alpha=alpha,
                                    keep = TRUE)
        
        lambda.best[outer.t] <- CV.inner.model$lambda.min
        beta[[outer.t]] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[-1]
        mu[outer.t] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[1]
        R2[outer.t] <- R2(mu[outer.t] + (cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]]) , Y.outer.test)
        RMSE[outer.t] <- RMSE(mu[outer.t] + (cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]]) , Y.outer.test)
        Rho[outer.t] <- cor(mu[outer.t] + (cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]]) , Y.outer.test)
        VarA[outer.t] <- var(cbind(X.outer.test,Z.outer.test)%*%beta[[outer.t]])
        VarE[outer.t] <- var(Y.outer.train)-VarA[outer.t]
        
      }
    }
    
    # Parallel Comp Output
    return(list( 
      "beta"= beta[[which.max(R2)]],
      "lambda" = lambda.best[which.max(R2)],
      "R2.list" = R2[which.max(R2)],
      "RMSE.list" = RMSE[which.max(R2)],
      "Rho.list" = Rho[which.max(R2)],
      "VarA.list" = VarA[which.max(R2)],
      "VarE.list" = VarE[which.max(R2)],
      "mu" = mu[which.max(R2)]
    ))
    
  }
  close(pb)
  stopCluster(cl)
  
  # Output
  return(list("beta" = do.call("rbind", Y_pred[which(1:length(Y_pred) %% 8 == 1)]), 
              "lambda" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 2)]), 
              "R2.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)]),
              "R2" = list("R2.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)])),
                          "R2.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)]))
              ),
              "RMSE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)]),
              "RMSE" = list("RMSE.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)])),
                            "RMSE.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)]))
              ),
              "Rho.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)]),
              "Rho" = list("Rho.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)])),
                           "Rho.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)]))
              ),
              "VarA.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 6)]),
              "VarE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 7)]),
              "mu" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 0)]),
              "aplha"=alpha))
}

# Nested Cross Validation Ridge Regression --------------------------------
#' @description Ridge Regression using the R package glmnet with standard nested cross validation framework
#' @param Y phenotype vector
#' @param First.Omic First omic matrix
#' @param Second.Omic Second omic matrix to concatenate with the First omic matrix
#' @param Third.Omic Third omic matrix to concatenate with the First and second omic matrix
#' @param cores cores number for parallel computing
#' @param alpha glmnet parameter alpha=0 for rigde regression
#' @param CV.fold Fold Sampling for Nested Cross Validation
#' @return list of : variable effects, mu, ridge reg lambda, 
#' @return model accuracies, variation du to genetic and variation du to errors
Ncv_RidgeReg_2 <- function(Y, First.Omic, Second.Omic=NULL, Third.Omic=NULL, cores = 1, alpha = 0, CV.fold){
  # Package Requirements
  require(doSNOW)
  require(glmnet)
  require(caret)
  
  # checks
  X <- First.Omic
  Z <- Second.Omic
  W <- Third.Omic
  stopifnot(length(Y) == nrow(X))
  if(!is.null(Z)){
    stopifnot(length(Y) == nrow(Z))
  }
  if(!is.null(W)){
    stopifnot(length(Y) == nrow(W))
  }
  sel <- !is.na(Y)
  X <- X[sel,]
  Y <- Y[sel]
  if(!is.null(Z)){
    Z <- Z[sel,]
  }
  if(!is.null(W)){
    W <- W[sel,]
  }
  
  # Nested cross-validation parameters
  iter <- length(CV.fold)
  outer.fold <- length(unique(CV.fold[[1]]$out.id))
  
  # Set parallel computing
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  # Computing Regression
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("glmnet","caret"), .options.snow=opts) %dopar% {
    
    # Output object
    beta <- list(); lambda.best <- mu <- R2 <- RMSE <- Rho <- VarA <- VarE <- numeric()
    
    # Outer-fold
    outer.foldid=CV.fold[[i]]$out.id
    
    # Outer-fold loop
    for(outer.t in 1:outer.fold){
      
      outer.test <- which(outer.foldid==outer.t)
      X.outer.test <- X[outer.test,]
      Y.outer.test <- Y[outer.test]
      outer.train <- which(outer.foldid!=outer.t)
      Y.outer.train <- Y[outer.train]
      X.outer.train <- X[outer.train,]
      if(!is.null(Z)){
        Z.outer.test <- Z[outer.test,]
        Z.outer.train <- Z[outer.train,]
      }
      if(!is.null(W)){
        W.outer.test <- W[outer.test,]
        W.outer.train <- W[outer.train,]
      }
      
      # Inner-fold
      inner.foldid=CV.fold[[i]]$in.id[[outer.t]]
      
      # Fit inner model using inner-fold CV: model
      if(is.null(Z) && is.null(W)){
        CV.inner.model <- cv.glmnet(x=X.outer.train, y=Y.outer.train, foldid=inner.foldid, alpha=alpha, keep = TRUE)
        
        lambda.best[outer.t] <- CV.inner.model$lambda.min
        beta[[outer.t]] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[-1]
        mu[outer.t] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[1]
        R2[outer.t] <- R2(mu[outer.t] + (X.outer.test %*% beta[[outer.t]]), Y.outer.test)
        RMSE[outer.t] <- RMSE(mu[outer.t] + (X.outer.test %*% beta[[outer.t]]), Y.outer.test)
        Rho[outer.t] <- cor(mu[outer.t] + (X.outer.test %*% beta[[outer.t]]), Y.outer.test)
        VarA[outer.t] <- var(X.outer.test %*% beta[[outer.t]])
        VarE[outer.t] <- var(Y.outer.train) - VarA[outer.t]
        
      } else if(!is.null(Z) && is.null(W)){
        CV.inner.model <- cv.glmnet(x=cbind(X.outer.train, Z.outer.train), y=Y.outer.train, foldid=inner.foldid, alpha=alpha, keep = TRUE)
        
        lambda.best[outer.t] <- CV.inner.model$lambda.min
        beta[[outer.t]] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[-1]
        mu[outer.t] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[1]
        R2[outer.t] <- R2(mu[outer.t] + (cbind(X.outer.test, Z.outer.test) %*% beta[[outer.t]]), Y.outer.test)
        RMSE[outer.t] <- RMSE(mu[outer.t] + (cbind(X.outer.test, Z.outer.test) %*% beta[[outer.t]]), Y.outer.test)
        Rho[outer.t] <- cor(mu[outer.t] + (cbind(X.outer.test, Z.outer.test) %*% beta[[outer.t]]), Y.outer.test)
        VarA[outer.t] <- var(cbind(X.outer.test, Z.outer.test) %*% beta[[outer.t]])
        VarE[outer.t] <- var(Y.outer.train) - VarA[outer.t]
        
      } else if(!is.null(Z) && !is.null(W)){
        CV.inner.model <- cv.glmnet(x=cbind(X.outer.train, Z.outer.train, W.outer.train), y=Y.outer.train, foldid=inner.foldid, alpha=alpha, keep = TRUE)
        
        lambda.best[outer.t] <- CV.inner.model$lambda.min
        beta[[outer.t]] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[-1]
        mu[outer.t] <- as.vector(coef(CV.inner.model, s = "lambda.min"))[1]
        R2[outer.t] <- R2(mu[outer.t] + (cbind(X.outer.test, Z.outer.test, W.outer.test) %*% beta[[outer.t]]), Y.outer.test)
        RMSE[outer.t] <- RMSE(mu[outer.t] + (cbind(X.outer.test, Z.outer.test, W.outer.test) %*% beta[[outer.t]]), Y.outer.test)
        Rho[outer.t] <- cor(mu[outer.t] + (cbind(X.outer.test, Z.outer.test, W.outer.test) %*% beta[[outer.t]]), Y.outer.test)
        VarA[outer.t] <- var(cbind(X.outer.test, Z.outer.test, W.outer.test) %*% beta[[outer.t]])
        VarE[outer.t] <- var(Y.outer.train) - VarA[outer.t]
        
      }
    }
    
    # Parallel Comp Output
  if(length(R2) == 0 || all(is.na(R2))) {
    print("R2 est vide ou contient uniquement des NA.")
    print(paste("Iteration:", i, "Outer fold:", outer.t))
    print(length(outer.foldid))
    print(length(inner.foldid))
    return(NULL)
  } else {
    return(list( 
      "beta" = beta[[which.max(R2)]],
      "lambda" = lambda.best[which.max(R2)],
      "R2.list" = R2[which.max(R2)],
      "RMSE.list" = RMSE[which.max(R2)],
      "Rho.list" = Rho[which.max(R2)],
      "VarA.list" = VarA[which.max(R2)],
      "VarE.list" = VarE[which.max(R2)],
      "mu" = mu[which.max(R2)]
    ))
   }
  }
  close(pb)
  stopCluster(cl)
  
  # Output
  return(list("beta" = do.call("rbind", Y_pred[which(1:length(Y_pred) %% 8 == 1)]), 
              "lambda" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 2)]), 
              "R2.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)]),
              "R2" = list("R2.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)])),
                          "R2.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 3)]))
              ),
              "RMSE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)]),
              "RMSE" = list("RMSE.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)])),
                            "RMSE.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 4)]))
              ),
              "Rho.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)]),
              "Rho" = list("Rho.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)])),
                           "Rho.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 5)]))
              ),
              "VarA.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 6)]),
              "VarE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 7)]),
              "mu" = do.call("c", Y_pred[which(1:length(Y_pred) %% 8 == 0)]),
              "alpha" = alpha))
}


Ncv_RFReg <- function(Y, trait, First.Omic, Second.Omic = NULL, Third.Omic = NULL, cores = 1, CV.fold, ntree = 500, mtry = NULL) {
    print_message("Random Forest", trait, 
                  First.Omic = deparse(substitute(First.Omic)),
                  Second.Omic = deparse(substitute(Second.Omic)),
                  Third.Omic = deparse(substitute(Third.Omic))
    )
  
  # Package Requirements
  require(doSNOW)
  require(randomForest)
  require(caret)
  
  # Checks
  X <- First.Omic
  Z <- Second.Omic
  W <- Third.Omic
  stopifnot(length(Y) == nrow(X))
  if (!is.null(Z)) {
    stopifnot(length(Y) == nrow(Z))
  }
  if (!is.null(W)) {
    stopifnot(length(Y) == nrow(W))
  }
  sel <- !is.na(Y)
  X <- X[sel, ]
  Y <- Y[sel]
  if (!is.null(Z)) {
    Z <- Z[sel, ]
  }
  if (!is.null(W)) {
    W <- W[sel, ]
  }
  
  # Nested cross-validation parameters
  iter <- length(CV.fold)
  outer.fold <- length(unique(CV.fold[[1]]$out.id))
  
  # Set parallel computing
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Computing Regression
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("randomForest", "caret"), .options.snow = opts) %dopar% {
    
    # Output object
    R2 <- RMSE <- Rho <- numeric()
    
    # Outer-fold
    outer.foldid = CV.fold[[i]]$out.id
    
    # Outer-fold loop
    for (outer.t in 1:outer.fold) {
      
      outer.test <- which(outer.foldid == outer.t)
      X.outer.test <- X[outer.test, ]
      Y.outer.test <- Y[outer.test]
      outer.train <- which(outer.foldid != outer.t)
      Y.outer.train <- Y[outer.train]
      X.outer.train <- X[outer.train, ]
      if (!is.null(Z)) {
        Z.outer.test <- Z[outer.test, ]
        Z.outer.train <- Z[outer.train, ]
      }
      if (!is.null(W)) {
        W.outer.test <- W[outer.test, ]
        W.outer.train <- W[outer.train, ]
      }
      
      # Combine omic data
      combined_train <- cbind(X.outer.train, Z.outer.train, if (!is.null(W)) W.outer.train else NULL)
      combined_test <- cbind(X.outer.test, Z.outer.test, if (!is.null(W)) W.outer.test else NULL)
      
      # Train Random Forest model
      rf_model <- randomForest(x = combined_train, y = Y.outer.train, ntree = ntree, mtry = ifelse(is.null(mtry), sqrt(ncol(combined_train)), mtry))
      
      # Predictions
      Y_pred_test <- predict(rf_model, combined_test)
      
      # Metrics
      R2[outer.t] <- R2(Y_pred_test, Y.outer.test)
      RMSE[outer.t] <- RMSE(Y_pred_test, Y.outer.test)
      Rho[outer.t] <- cor(Y_pred_test, Y.outer.test)
    }
    
    # Parallel Comp Output
    return(list(
      "R2.list" = R2[which.max(R2)],
      "RMSE.list" = RMSE[which.max(R2)],
      "Rho.list" = Rho[which.max(R2)]
    ))
    
  }
  close(pb)
  stopCluster(cl)
  
  # Output
  return(list(
    "R2.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)]),
    "R2" = list(
      "R2.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)])),
      "R2.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)]))
    ),
    "RMSE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)]),
    "RMSE" = list(
      "RMSE.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)])),
      "RMSE.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)]))
    ),
    "Rho.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)]),
    "Rho" = list(
      "Rho.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)])),
      "Rho.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)]))
    )
  ))
}



Ncv_CNNReg <- function(Y, trait, First.Omic, Second.Omic = NULL, Third.Omic = NULL, cores = 1, CV.fold, 
                       epochs = 50, batch_size = 32, learning_rate = 0.001) {
    print_message("CNN", trait, 
                  First.Omic = deparse(substitute(First.Omic)),
                  Second.Omic = deparse(substitute(Second.Omic)),
                  Third.Omic = deparse(substitute(Third.Omic))
    )
  
  # Package Requirements
  require(doSNOW)
  require(keras)
  require(caret)
  
  # Checks
  X <- First.Omic
  Z <- Second.Omic
  W <- Third.Omic
  stopifnot(length(Y) == nrow(X))
  if (!is.null(Z)) {
    stopifnot(length(Y) == nrow(Z))
  }
  if (!is.null(W)) {
    stopifnot(length(Y) == nrow(W))
  }
  sel <- !is.na(Y)
  X <- X[sel, ]
  Y <- Y[sel]
  if (!is.null(Z)) {
    Z <- Z[sel, ]
  }
  if (!is.null(W)) {
    W <- W[sel, ]
  }
  
  # Nested cross-validation parameters
  iter <- length(CV.fold)
  outer.fold <- length(unique(CV.fold[[1]]$out.id))
  
  # Set parallel computing
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Function to build CNN model
  build_cnn_model <- function(input_shape) {
    model <- keras_model_sequential() %>%
      layer_conv_1d(filters = 32, kernel_size = 3, activation = 'relu', input_shape = input_shape) %>%
      layer_max_pooling_1d(pool_size = 2) %>%
      layer_conv_1d(filters = 64, kernel_size = 3, activation = 'relu') %>%
      layer_max_pooling_1d(pool_size = 2) %>%
      layer_flatten() %>%
      layer_dense(units = 64, activation = 'relu') %>%
      layer_dense(units = 1)
    
    model %>% compile(
      optimizer = optimizer_adam(learning_rate = learning_rate),
      loss = 'mse',
      metrics = c('mean_squared_error')
    )
    
    return(model)
  }
  
  # Computing Regression
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("keras", "caret"), .options.snow = opts) %dopar% {
    
    # Output object
    R2 <- RMSE <- Rho <- numeric()
    
    # Outer-fold
    outer.foldid = CV.fold[[i]]$out.id
    
    # Outer-fold loop
    for (outer.t in 1:outer.fold) {
      
      outer.test <- which(outer.foldid == outer.t)
      X.outer.test <- X[outer.test, ]
      Y.outer.test <- Y[outer.test]
      outer.train <- which(outer.foldid != outer.t)
      Y.outer.train <- Y[outer.train]
      X.outer.train <- X[outer.train, ]
      if (!is.null(Z)) {
        Z.outer.test <- Z[outer.test, ]
        Z.outer.train <- Z[outer.train, ]
      }
      if (!is.null(W)) {
        W.outer.test <- W[outer.test, ]
        W.outer.train <- W[outer.train, ]
      }
      
      # Combine omic data
      combined_train <- cbind(X.outer.train, Z.outer.train, if (!is.null(W)) W.outer.train else NULL)
      combined_test <- cbind(X.outer.test, Z.outer.test, if (!is.null(W)) W.outer.test else NULL)
      
      # Reshape for CNN (samples, time_steps, features)
      combined_train <- array_reshape(combined_train, c(nrow(combined_train), ncol(combined_train), 1))
      combined_test <- array_reshape(combined_test, c(nrow(combined_test), ncol(combined_test), 1))
      
      # Build and train CNN model
      cnn_model <- build_cnn_model(input_shape = c(ncol(combined_train), 1))
      cnn_model %>% fit(
        x = combined_train, y = Y.outer.train,
        epochs = epochs, batch_size = batch_size, verbose = 0
      )
      
      # Predictions
      Y_pred_test <- predict(cnn_model, combined_test)
      
      # Metrics
      R2[outer.t] <- R2(as.vector(Y_pred_test), Y.outer.test)
      RMSE[outer.t] <- RMSE(as.vector(Y_pred_test), Y.outer.test)
      Rho[outer.t] <- cor(as.vector(Y_pred_test), Y.outer.test)
    }
    
    # Parallel Comp Output
    return(list(
      "R2.list" = R2[which.max(R2)],
      "RMSE.list" = RMSE[which.max(R2)],
      "Rho.list" = Rho[which.max(R2)]
    ))
    
  }
  close(pb)
  stopCluster(cl)
  
  # Output
  return(list(
    "R2.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)]),
    "R2" = list(
      "R2.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)])),
      "R2.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)]))
    ),
    "RMSE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)]),
    "RMSE" = list(
      "RMSE.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)])),
      "RMSE.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)]))
    ),
    "Rho.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)]),
    "Rho" = list(
      "Rho.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)])),
      "Rho.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)]))
    )
  ))
}



Ncv_AlexNETReg <- function(Y, trait,First.Omic, Second.Omic = NULL, Third.Omic = NULL, cores = 1, CV.fold, 
                       epochs = 50, batch_size = 32, learning_rate = 0.001) {

  print_message("AlexNET", trait, 
                First.Omic = deparse(substitute(First.Omic)),
                Second.Omic = deparse(substitute(Second.Omic)),
                Third.Omic = deparse(substitute(Third.Omic))
                )  # Package Requirements
  require(doSNOW)
  require(keras)
  require(caret)
  
  # Checks
  X <- First.Omic
  Z <- Second.Omic
  W <- Third.Omic
  stopifnot(length(Y) == nrow(X))
  if (!is.null(Z)) {
    stopifnot(length(Y) == nrow(Z))
  }
  if (!is.null(W)) {
    stopifnot(length(Y) == nrow(W))
  }
  sel <- !is.na(Y)
  X <- X[sel, ]
  Y <- Y[sel]
  if (!is.null(Z)) {
    Z <- Z[sel, ]
  }
  if (!is.null(W)) {
    W <- W[sel, ]
  }
  
  # Nested cross-validation parameters
  iter <- length(CV.fold)
  outer.fold <- length(unique(CV.fold[[1]]$out.id))
  
  # Set parallel computing
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Function to build AlexNet model
  build_alexnet_model <- function(input_shape) {
    model <- keras_model_sequential() %>%
      layer_conv_2d(filters = 96, kernel_size = c(11, 11), strides = c(4, 4), activation = 'relu', 
                    input_shape = input_shape) %>%
      layer_max_pooling_2d(pool_size = c(3, 3), strides = c(2, 2)) %>%
      layer_conv_2d(filters = 256, kernel_size = c(5, 5), padding = 'same', activation = 'relu') %>%
      layer_max_pooling_2d(pool_size = c(3, 3), strides = c(2, 2)) %>%
      layer_conv_2d(filters = 384, kernel_size = c(3, 3), padding = 'same', activation = 'relu') %>%
      layer_conv_2d(filters = 384, kernel_size = c(3, 3), padding = 'same', activation = 'relu') %>%
      layer_conv_2d(filters = 256, kernel_size = c(3, 3), padding = 'same', activation = 'relu') %>%
      layer_max_pooling_2d(pool_size = c(3, 3), strides = c(2, 2)) %>%
      layer_flatten() %>%
      layer_dense(units = 4096, activation = 'relu') %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 4096, activation = 'relu') %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 1)
    
    model %>% compile(
      optimizer = optimizer_adam(learning_rate = learning_rate),
      loss = 'mse',
      metrics = c('mean_squared_error')
    )
    
    return(model)
  }
  
  # Computing Regression
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("keras", "caret"), .options.snow = opts) %dopar% {
    
    # Output object
    R2 <- RMSE <- Rho <- numeric()
    
    # Outer-fold
    outer.foldid = CV.fold[[i]]$out.id
    
    # Outer-fold loop
    for (outer.t in 1:outer.fold) {
      
      outer.test <- which(outer.foldid == outer.t)
      X.outer.test <- X[outer.test, ]
      Y.outer.test <- Y[outer.test]
      outer.train <- which(outer.foldid != outer.t)
      Y.outer.train <- Y[outer.train]
      X.outer.train <- X[outer.train, ]
      if (!is.null(Z)) {
        Z.outer.test <- Z[outer.test, ]
        Z.outer.train <- Z[outer.train, ]
      }
      if (!is.null(W)) {
        W.outer.test <- W[outer.test, ]
        W.outer.train <- W[outer.train, ]
      }
      
      # Combine omic data
      combined_train <- cbind(X.outer.train, Z.outer.train, if (!is.null(W)) W.outer.train else NULL)
      combined_test <- cbind(X.outer.test, Z.outer.test, if (!is.null(W)) W.outer.test else NULL)
      
      # Reshape for AlexNet (samples, height, width, channels)
      combined_train <- array_reshape(combined_train, c(nrow(combined_train), ncol(combined_train), 1, 1))
      combined_test <- array_reshape(combined_test, c(nrow(combined_test), ncol(combined_test), 1, 1))
      
      # Build and train AlexNet model
      alexnet_model <- build_alexnet_model(input_shape = c(ncol(combined_train), 1, 1))
      alexnet_model %>% fit(
        x = combined_train, y = Y.outer.train,
        epochs = epochs, batch_size = batch_size, verbose = 0
      )
      
      # Predictions
      Y_pred_test <- predict(alexnet_model, combined_test)
      
      # Metrics
      R2[outer.t] <- R2(as.vector(Y_pred_test), Y.outer.test)
      RMSE[outer.t] <- RMSE(as.vector(Y_pred_test), Y.outer.test)
      Rho[outer.t] <- cor(as.vector(Y_pred_test), Y.outer.test)
    }
    
    # Parallel Comp Output
    return(list(
      "R2.list" = R2[which.max(R2)],
      "RMSE.list" = RMSE[which.max(R2)],
      "Rho.list" = Rho[which.max(R2)]
    ))
    
  }
  close(pb)
  stopCluster(cl)
  
  # Output
  return(list(
    "R2.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)]),
    "R2" = list(
      "R2.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)])),
      "R2.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 1)]))
    ),
    "RMSE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)]),
    "RMSE" = list(
      "RMSE.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)])),
      "RMSE.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 2)]))
    ),
    "Rho.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)]),
    "Rho" = list(
      "Rho.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)])),
      "Rho.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 3 == 0)]))
    )
  ))
}


Ncv_AlexNETReg_up <- function(Y, First.Omic, Second.Omic = NULL, Third.Omic = NULL, cores = 1, CV.fold, 
                       epochs = 50, batch_size = 32, learning_rate = 0.001) {
  # Package Requirements
  require(doSNOW)
  require(keras)
  require(caret)
  
  # Checks
  X <- First.Omic
  Z <- Second.Omic
  W <- Third.Omic
  stopifnot(length(Y) == nrow(X))
  if (!is.null(Z)) {
    stopifnot(length(Y) == nrow(Z))
  }
  if (!is.null(W)) {
    stopifnot(length(Y) == nrow(W))
  }
  sel <- !is.na(Y)
  X <- X[sel, ]
  Y <- Y[sel]
  if (!is.null(Z)) {
    Z <- Z[sel, ]
  }
  if (!is.null(W)) {
    W <- W[sel, ]
  }
  
  # Nested cross-validation parameters
  iter <- length(CV.fold)
  outer.fold <- length(unique(CV.fold[[1]]$out.id))
  
  # Set parallel computing
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Function to build AlexNet model
  build_alexnet_model <- function(input_shape) {
    model <- keras_model_sequential() %>%
      layer_conv_2d(filters = 96, kernel_size = c(11, 11), strides = c(4, 4), activation = 'relu', 
                    input_shape = input_shape) %>%
      layer_max_pooling_2d(pool_size = c(3, 3), strides = c(2, 2)) %>%
      layer_conv_2d(filters = 256, kernel_size = c(5, 5), padding = 'same', activation = 'relu') %>%
      layer_max_pooling_2d(pool_size = c(3, 3), strides = c(2, 2)) %>%
      layer_conv_2d(filters = 384, kernel_size = c(3, 3), padding = 'same', activation = 'relu') %>%
      layer_conv_2d(filters = 384, kernel_size = c(3, 3), padding = 'same', activation = 'relu') %>%
      layer_conv_2d(filters = 256, kernel_size = c(3, 3), padding = 'same', activation = 'relu') %>%
      layer_max_pooling_2d(pool_size = c(3, 3), strides = c(2, 2)) %>%
      layer_flatten() %>%
      layer_dense(units = 4096, activation = 'relu') %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 4096, activation = 'relu') %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 1)
    
    model %>% compile(
      optimizer = optimizer_adam(learning_rate = learning_rate),
      loss = 'mse',
      metrics = c('mean_squared_error')
    )
    
    return(model)
  }
  
  # Function to calculate Grad-CAM
  calculate_gradcam <- function(model, input_data, target_index = NULL) {
    grad_model <- keras_model(inputs = model$input, outputs = list(model$get_layer(index = length(model$layers) - 1)$output, 
                                                                   model$output))
    
    with(tf$GradientTape() %as% tape, {
      inputs <- tf$convert_to_tensor(input_data)
      tape$watch(inputs)
      outputs <- grad_model(inputs)
      conv_outputs <- outputs[[1]]
      predictions <- outputs[[2]]
      loss <- predictions[, target_index]
    })
    
    grads <- tape$gradient(loss, conv_outputs)
    pooled_grads <- tf$reduce_mean(grads, axis = c(1, 2))
    conv_outputs <- conv_outputs[1,,,]
    pooled_grads <- pooled_grads[1,]
    
    for (i in seq_along(pooled_grads)) {
      conv_outputs[,,i] <- conv_outputs[,,i] * pooled_grads[i]
    }
    
    heatmap <- tf$reduce_mean(conv_outputs, axis = -1)$numpy()
    heatmap <- pmax(heatmap, 0)
    heatmap <- heatmap / max(heatmap)
    return(heatmap)
  }
  
  # Computing Regression
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("keras", "caret"), .options.snow = opts) %dopar% {
    
    # Output object
    R2 <- RMSE <- Rho <- list()
    GradCAM_outputs <- list()
    
    # Outer-fold
    outer.foldid = CV.fold[[i]]$out.id
    
    # Outer-fold loop
    for (outer.t in 1:outer.fold) {
      
      outer.test <- which(outer.foldid == outer.t)
      X.outer.test <- X[outer.test, ]
      Y.outer.test <- Y[outer.test]
      outer.train <- which(outer.foldid != outer.t)
      Y.outer.train <- Y[outer.train]
      X.outer.train <- X[outer.train, ]
      if (!is.null(Z)) {
        Z.outer.test <- Z[outer.test, ]
        Z.outer.train <- Z[outer.train, ]
      }
      if (!is.null(W)) {
        W.outer.test <- W[outer.test, ]
        W.outer.train <- W[outer.train, ]
      }
      
      # Combine omic data
      combined_train <- cbind(X.outer.train, Z.outer.train, if (!is.null(W)) W.outer.train else NULL)
      combined_test <- cbind(X.outer.test, Z.outer.test, if (!is.null(W)) W.outer.test else NULL)
      
      # Reshape for AlexNet (samples, height, width, channels)
      combined_train <- array_reshape(combined_train, c(nrow(combined_train), ncol(combined_train), 1, 1))
      combined_test <- array_reshape(combined_test, c(nrow(combined_test), ncol(combined_test), 1, 1))
      
      # Build and train AlexNet model
      alexnet_model <- build_alexnet_model(input_shape = c(ncol(combined_train), 1, 1))
      alexnet_model %>% fit(
        x = combined_train, y = Y.outer.train,
        epochs = epochs, batch_size = batch_size, verbose = 0
      )
      
      # Predictions
      Y_pred_test <- predict(alexnet_model, combined_test)
      
      # Metrics
      R2[[outer.t]] <- R2(as.vector(Y_pred_test), Y.outer.test)
      RMSE[[outer.t]] <- RMSE(as.vector(Y_pred_test), Y.outer.test)
      Rho[[outer.t]] <- cor(as.vector(Y_pred_test), Y.outer.test)
      
      # Grad-CAM Calculation
      GradCAM_outputs[[outer.t]] <- calculate_gradcam(alexnet_model, combined_test, target_index = 1)
    }
    
    # Parallel Comp Output
    return(list(
      "R2.list" = R2[which.max(R2)],
      "RMSE.list" = RMSE[which.max(R2)],
      "Rho.list" = Rho[which.max(R2)],
      "GradCAM" = GradCAM_outputs
    ))
    
  }
  close(pb)
  stopCluster(cl)
  
  # Output
  return(list(
    "R2.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 1)]),
    "R2" = list(
      "R2.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 1)])),
      "R2.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 1)]))
    ),
    "RMSE.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 2)]),
    "RMSE" = list(
      "RMSE.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 2)])),
      "RMSE.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 2)]))
    ),
    "Rho.list" = do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 3)]),
    "Rho" = list(
      "Rho.mean" = mean(do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 3)])),
      "Rho.sd" = sd(do.call("c", Y_pred[which(1:length(Y_pred) %% 4 == 3)]))
    ),
    "GradCAM" = Y_pred[which(1:length(Y_pred) %% 4 == 0)]
  ))
}
