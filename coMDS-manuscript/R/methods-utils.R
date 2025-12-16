#' Ensemble Distances in Meta-Sepc Visualization
#'
#' @description Copied (and renamed) from `ensemble.viz()` located at
#'   https://github.com/rongstat/meta-visualization/blob/main/R%20Codes/main_fun.R
#'
#' @param data_list A list of 2-dimensional embeddings
#' @param method_names Names of the candidate visualizations
#' @param original_data An option to use the original data for the quality
#'   assessment, instead of using eigenscores
#'
#' @returns A list containing:
#' - ensemble_dist_mat: a meta-distance for meta-visualization
#' - eigenscore: eigenscores for candidate visualizations
#' - method_names: names of the method ordered by averaged eigenscores.
ensemble_metaspec <- function(data_list, method_names = NA, original_data = NA,
                              verbose = FALSE) {
  n <- dim(data_list[[1]])[1]
  K <- length(data_list)
  if (is.na(method_names)) {
    method_names <- names(data_list)
  }

  if (is.na(original_data)) {
    ########## obtain weights
    ensemble_mat <- matrix(ncol = n, nrow = n)
    weight <- matrix(ncol = K, nrow = n)
    for (j in 1:n) {
      local_dist <- matrix(ncol = n, nrow = K)
      for (i in 1:K) {
        local_dist[i, ] <- sqrt(
          rowSums(t(data_list[[i]][j, ] - t(data_list[[i]]))^2)
        )
      }
      comp_mat <- matrix(ncol = K, nrow = K)
      embed_mat_norm <- list()
      for (i in 1:K) {
        for (k in 1:K) {
          comp_mat[i, k] <- sum(local_dist[i, ] * local_dist[k, ]) /
            sqrt(sum(local_dist[k, ]^2)) / sqrt(sum(local_dist[i, ]^2))
        }
        embed_mat_norm[[i]] <- local_dist[i, ] / sqrt(sum(local_dist[i, ]^2))
      }
      weight[j, ] <- abs(eigen(comp_mat)$vectors[, 1])
      ensemble_mat[, j] <- apply(
        do.call(cbind, embed_mat_norm), 1, weighted.mean, w = weight[j, ]
      ) * sum(weight[j, ])
      if (verbose) {
        if (j / 1000 == floor(j / 1000)) {
          print(paste0(j, " samples done!"))
        }
      }
    }
  } else {
    ensemble_mat <- matrix(ncol = n, nrow = n)
    weight <- matrix(ncol = K, nrow = n)
    data_mat <- as.matrix(dist(original_data))
    for (j in 1:n) {
      local_dist <- matrix(ncol = n, nrow = K)
      for (i in 1:K) {
        local_dist[i, ] <- sqrt(
          rowSums(t(data_list[[i]][j, ] - t(data_list[[i]]))^2)
        )
      }
      eigen_score <- c()
      embed_mat_norm <- list()
      for (i in 1:K) {
        eigen_score[i] <- sum(data_mat[, j] * local_dist[i, ]) /
          sqrt(sum(data_mat[, j]^2)) / sqrt(sum(local_dist[i, ]^2))
        embed_mat_norm[[i]] <- local_dist[i, ] / sqrt(sum(local_dist[i, ]^2))
      }
      eigen_score[which(eigen_score < 0)] <- 0
      weight[j, ] <- eigen_score^2 / sum(eigen_score^2)
      ensemble_mat[, j] <- apply(
        do.call(cbind, embed_mat_norm), 1, weighted.mean, w = weight[j, ]
      ) * sum(weight[j, ])
    }
  }
  return(list(
    ensemble_dist_mat = (ensemble_mat + t(ensemble_mat)) / 2,
    eigenscore = weight,
    method_names = method_names[order(colMeans(weight), decreasing = T)]
  ))
}
