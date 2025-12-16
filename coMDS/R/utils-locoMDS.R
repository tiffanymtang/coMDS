#' Get neighborhood indicators for LoCoMDS
#' @keywords internal
get_neighbors_list <- function(delta, percentile) {
  M <- length(delta)
  dist_vecs <- list()
  for (m in 1:M) {
    dist_vecs[[m]] <- c(as.dist(delta[[m]]))
  }
  eps_vec <- unlist(lapply(delta, quantile, probs = percentile), use.names = FALSE)
  result <- get_neighbors(dist_vecs, eps_vec)
  return(result)
}


#' Get weight matrices based upon local neighborhoods
#' @keywords internal
init_local_weights <- function(neighbor_list) {
  M <- length(neighbor_list)
  N <- length(neighbor_list[[1]])
  weights <- repList(matrix(1e-12, N, N), M)
  for (m in 1:M) {
    this_weight <- weights[[m]]
    this_neighbor_list <- neighbor_list[[m]]
    for (i in 1:N) {
      this_neighbors <- this_neighbor_list[[i]]
      for (j in this_neighbors) {
        j <- j + 1
        this_weight[i, j] <- 1
      }
    }
    weights[[m]] <- this_weight
  }
  return(weights)
}


#' Get t value for regularization in LoCoMDS
#' @keywords internal
get_t <- function(delta, weights, tau, percentile) {
  dist_vecs <- purrr::map2(delta, weights, ~ c(.x[.y == 1]))
  coeff <- ifelse(percentile == 1, 1, (percentile) / (1 - percentile))
  result <- coeff * median(unlist(dist_vecs)) * tau
  return(result)
}

#' Get t values for regularization in LoCoMDS
#' @keywords internal
get_ts <- function(delta, weights, tau, percentile) {
  dist_vecs <- purrr::map2(delta, weights, ~ c(.x[.y == 1]))
  coeff <- ifelse(percentile == 1, 1, (percentile) / (1 - percentile))
  result <- coeff * sapply(dist_vecs, median) * tau
  return(result)
}


#' Get B matrix for LoCoMDS
#' @keywords internal
loco_bmat <- function(diss, wgths, d, tval, ones_matrix, eps = 1E-12) {
  z <- ifelse(d < eps, 1, 0)
  b <- as.matrix((wgths * diss + 0.5 * tval * (ones_matrix - wgths)) * (1 - z) / (d + z))
  r <- rowSums(b)
  return(diag(r) - b)
}


#' Plot LCMC
#' @keywords internal
plot_locoMDS_lcmc <- function(lcmc_result) {
  plt <- lcmc_result |>
    dplyr::mutate(
      parameter_name = sprintf("tau = %s | percentile = %s", tau, percentile),
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = k,
      y = adjusted_lcmc,
      shape = parameter_name,
      color = parameter_name
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Neighborhood Size (k)", y = "Adjusted LCMC",
      color = "Hyperparameters", shape = "Hyperparameters"
    ) +
    ggplot2::theme_minimal()
  return(plt)
}
