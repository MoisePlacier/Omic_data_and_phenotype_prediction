#' rank_features
#'
#' @description ranks features based on their importance in predicting a response variable using Sparse Partial Least Squares (sPLS).
#'
#' @param X A matrix or data frame of predictor variables (features).
#' @param Y A vector or matrix of response variables.
#' @param ncomp The number of components to use in the sPLS model.
#' @param n_features The number of features to select based on their importance.
#'
#' @return A data frame of selected features ranked by their importance.
rank_features <- function(X, Y, ncomp = 2, n_features = 50) {
  X <- scale(X)
  Y <- scale(Y)
  if (n_features > 100) {
    n_features <- 100
  }

  res <- spls(X, Y, ncomp = ncomp, keepX = c(n_features, n_features))

  selected_genes <- selectVar(res, comp = 1)$X$value
  return(selected_genes)
}

#' vis_feature_ranking
#'
#' @description Visualizes the feature ranking results.
#'
#' @param X A matrix or data frame of predictor variables (features).
#' @param Y A vector or matrix of response variables.
#' @param ncomp The number of components to use in the sPLS model.
#' @param n_features The number of features to select based on their importance.
#'
#' @return A list containing the selected feature names.
select_features <- function(X, Y, ncomp = 2, n_features = 50) {
  X <- scale(X)
  Y <- scale(Y)
  max_features <- 100
  n_iter <- n_features %/% max_features
  remainder <- n_features %% max_features
  selected_genes <- list()
  if (n_iter) {
    for (i in 1:n_iter) {
      res <- spls(X, Y, ncomp = ncomp, keepX = c(max_features, max_features))
      c <- selectVar(res, comp = 1)$X$value
      selected_genes_i <- selectVar(res, comp = 1)$X$value
      selected_genes[[i]] <- selected_genes_i
      X <- X[, !colnames(X) %in% rownames(selected_genes_i)]
    }
  }

  if (remainder) {
    res <- spls(X, Y, ncomp = ncomp, keepX = c(remainder, remainder))
    selected_genes[[n_iter + 1]] <- selectVar(res, comp = 1)$X$value
  }

  # merge all selected genes into a single data frame
  selected_genes <- do.call(rbind, selected_genes)
  selected_genes <- as.data.frame(selected_genes)

  return(rownames(selected_genes))
}

#' feature_rank_per_omic
#'
#' @description Ranks features for each trait in the PoplarData object using the rank_features function.
#'
#' @param PoplarData A list containing omic data and phenotypes.
#' @param omic_data The specific omic data to process.
#' @param verbose An integer indicating the verbosity level of the output.
#'
#' @return A list of selected genes for each trait.
feature_rank_per_omic <- function(PoplarData, omic_data, ncomp = 2, n_features = 50) {
  res <- list()
  i <- 0
  for (trait in colnames(PoplarData$Phenotype)) {
    i <- i + 1
    selected_genes <- rank_features(PoplarData$Transcriptomic,
      PoplarData$Phenotype[[trait]],
      ncomp = ncomp,
      n_features = n_features
    )
    res[[trait]] <- selected_genes
  }
  return(res)
}

#' feature_ranking
#'
#' @description Ranks features for all omic data in the PoplarData object.
#'
#' @param PoplarData A list containing omic data and phenotypes.
#' @param verbose An integer indicating the verbosity level of the output.
#'
#' @return A list of feature rankings for each omic data type in the PoplarData object.
feature_ranking <- function(PoplarData, ncomp = 2, n_features = 50) {
  PoplarData.feature.ranking <- list()
  for (omic_data in names(PoplarData)) {
    if (omic_data == "Phenotype") nex
    feature_rank <- feature_rank_per_omic(PoplarData,
      PoplarData[[omic_data]],
      ncomp = ncomp,
      n_features = n_features
    )
    PoplarData.feature.ranking[[omic_data]] <- feature_rank
  }
  return(PoplarData.feature.ranking)
}

#' vis_selected_genes
#'
#' @description Visualizes the selected genes from the feature ranking.
#'
#' @param selected_genes A data frame of selected genes with their importance values.
#'
#' @return A ggplot object visualizing the selected genes.
vis_selected_genes <- function(selected_genes) {
  library(ggplot2)
  library(dplyr)

  selected_genes$gene <- rownames(selected_genes)
  selected_genes <- selected_genes %>%
    mutate(sign = ifelse(value.var >= 0, "Positive", "Negative"))

  selected_genes$value.var <- abs(selected_genes$value.var)

  selected_genes <- selected_genes %>%
    arrange(desc(abs(value.var))) %>%
    mutate(gene = factor(gene, levels = gene))

  ggplot(selected_genes, aes(x = gene, y = value.var, fill = sign)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +
    labs(x = "Gene", y = "Value", fill = "Sign") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
