#' Evaluation wrapper function
evaluate_dr_methods <- function(embed_list, data, labels,
                                metrics = c("triplet", "spearman", "rf", "lcmc_local"),
                                integration_method_start_index = 17,
                                n_replicates = 100) {
  metrics <- match.arg(metrics, several.ok = TRUE)
  n <- nrow(data)

  # get original distance matrix and new distance matrices
  ori_distance <- as.matrix(dist(data))
  embed_distance_ls <- purrr::map(embed_list, ~ as.matrix(dist(.x)))

  # evaluate metrics for each method
  eval_out <- purrr::map2(
    embed_list, embed_distance_ls,
    function(embed_data, embed_dist) {
      metrics_out <- list()
      for (metric in metrics) {
        if (metric == "triplet") {
          results <- replicate(
            n_replicates,
            { eval_triplet(ori_distance, embed_dist) }
          )
          metrics_out[[sprintf("%s_raw", metric)]] <- list(results)
          metrics_out[[metric]] <- mean(results)
        } else if (metric == "spearman") {
          results <- replicate(
            n_replicates,
            { eval_spearman(ori_distance, embed_dist) }
          )
          metrics_out[[sprintf("%s_raw", metric)]] <- list(results)
          metrics_out[[metric]] <- mean(results)
        } else if (metric == "lcmc_local") {
          local_ks <- seq(2, 20, 3)
          results <- coMDS::lcmc(data, embed_data, ks = local_ks)
          metrics_out[[sprintf("%s_raw", metric)]] <- list(results)
          metrics_out[[metric]] <- mean(results$lcmc)
        }
      }
      return(tibble::as_tibble(metrics_out))
    }
  ) |>
  dplyr::bind_rows(.id = ".method_name")
  if ("rf" %in% metrics) {
    rf_results <- purrr::map(
      1:(n_replicates*5), 
      ~ suppressMessages(eval_rf(embed_list, labels, integration_method_start_index))
    )
    rf_results <- purrr::map(rf_results, ~ purrr::map(.x, as.data.frame))
    rf_list <- purrr::map(1:length(rf_results[[1]]), function(j) {
      purrr::list_rbind(purrr::map(rf_results, ~ .x[[j]]))
    })
    eval_out <- eval_out |>
    dplyr::mutate(
      rf_raw = rf_list,
      rf     = purrr::map(rf_list, ~ as.data.frame(t(colMeans(.x, na.rm = TRUE))))
    )
  }
  return(eval_out)
}

#' Evaluate Random Triplet Accuracy
#'
#' @param ori_distance distance matrix for raw data/embeddings in original space.
#' @param new_distance distance matrix for raw data/embeddings in new (dimension-reduced) space.
#' @param num_triplets number of random triplet comparisons per sample. Defaults to 5.
#'
#' @return A numeric scalar: the mean proportion of triplet relationships preserved across all samples.
eval_triplet <- function(ori_distance, new_distance, num_triplets = 5) {
  n_samples <- nrow(ori_distance)
  result <- c()
  for (i in 1:n_samples) {
    this_accuracy <- 0
    triplet_index <- replicate(num_triplets, sample(c(1:n_samples)[-i], 2))
    for (j in 1:num_triplets) {
      this_j1 <- triplet_index[1, j]
      this_j2 <- triplet_index[2, j]
      ori_flag <- (ori_distance[i, this_j1] < ori_distance[i, this_j2])
      new_flag <- (new_distance[i, this_j1] < new_distance[i, this_j2])
      this_accuracy <- this_accuracy + (ori_flag == new_flag)
    }
    this_accuracy <- this_accuracy / num_triplets
    result <- c(result, this_accuracy)
  }
  result <- mean(result)
  return(result)
}


#' Evaluate Spearman Correlation between Two Distance Matrices
#'
#' @param ori_distance distance matrix for raw data/embeddings in original space.
#' @param new_distance distance matrix for raw data/embeddings in new (dimension-reduced) space.
#' @param percentage A numeric value (0–1) indicating what percentage of samples to use in the comparison. Defaults to 0.6.
#' @param max_n Maximum number of points to sample, even if the percentage exceeds this. Defaults to 800.
#'
#' @return A numeric scalar: the mean spearman correlation between the sample distances in original space and that in new space.
eval_spearman <- function(ori_distance, new_distance, percentage = 0.6, max_n = 800) {
  n_samples <- nrow(ori_distance)
  n_points <- min(floor(percentage * n_samples), max_n)
  sample_index <- sample(1:n_samples, n_points, replace = FALSE)
  dist_high <- as.vector(as.dist(ori_distance[sample_index, sample_index]))
  dist_low <- as.vector(as.dist(new_distance[sample_index, sample_index]))
  result <- cor(dist_high, dist_low, method = "spearman")
  return(result)
}


#' Evaluate Random Forest Classifier
#'
#' @param embed_list embeddings list.
#' @param labels true cluster labels.
#' @param integration_method_start_index index of the first integration method.
#' @param n_splits Number of folds for train test split. Defaults to 5.
#'
#' @return A list containing:
#' - `pred_accuracy`: prediction accuracy of the random forest classifier.
#' - `match_prob`: average probability of correct predictions.
#' - `diff_prob`: average probability of misclassifications.
eval_rf <- function(embed_list, labels, integration_method_start_index, n_splits = 5) {

  minMaxNorm <- function(x) (x - min(x)) / (max(x) - min(x))

  train_ind <- caret::createFolds(labels, k = n_splits, list = TRUE, returnTrain = TRUE)[[1]]

  rf_out <- lapply(embed_list, function(emb) {
    normed_coor <- data.frame(
      x1 = minMaxNorm(emb[, 1]),
      x2 = minMaxNorm(emb[, 2]),
      cluster = as.factor(labels)
    )

    train <- normed_coor[train_ind, ]
    test <- normed_coor[-train_ind, ]

    rf <- randomForest::randomForest(cluster ~ x1 + x2, data = train, ntree = 300)
    prediction_labels <- predict(rf, newdata = test)
    prediction_prob <- predict(rf, newdata = test, type = "prob")

    list(test = test, prediction_labels = prediction_labels, prediction_prob = prediction_prob)
  })

  selected_rf <- rf_out[integration_method_start_index:length(rf_out)]
  misclassified_pts_list <- lapply(selected_rf, function(x) {
    which(x$test$cluster != x$prediction_labels)
  })
  shared_misclassified <- Reduce(intersect, misclassified_pts_list)

  results <- lapply(rf_out, function(x) {
    test <- x$test
    pred <- x$prediction_labels
    prob <- x$prediction_prob

    # Accuracy
    pred_accuracy <- mean(test$cluster == pred)

    # Match probability
    match_idx <- which(test$cluster == pred)
    match_prob <- mean(mapply(function(r, c) prob[r, c], match_idx, pred[match_idx]), na.rm = TRUE)

    # Misclassification probability (shared)
    diff_idx <- shared_misclassified
    diff_prob <- mean(mapply(function(r, c) prob[r, c], diff_idx, pred[diff_idx]), na.rm = TRUE)

    # AUROC
    aucs <- sapply(levels(test$cluster), function(cls) {
      actual <- as.numeric(test$cluster == cls)
      preds <- prob[, cls]
      if (length(unique(actual)) < 2) NA else pROC::roc(actual, preds)$auc
    })
    macro_roc <- mean(aucs, na.rm = TRUE)
    multi_roc <- as.numeric(pROC::multiclass.roc(test$cluster, prob)$auc)

    list(
      pred_accuracy = pred_accuracy,
      match_prob = match_prob,
      diff_prob = diff_prob,
      macro_roc = macro_roc,
      multi_roc = multi_roc
    )
  })

  return(results)
}
