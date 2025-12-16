#' Evaluate Local Continuity Meta-Criterion (LCMC)
#'
#' @param x Raw data/embeddings in original space.
#' @param y Raw data/embeddings in new (dimension-reduced) space.
#' @param ks Vector of neighborhood sizes to evaluate LCMC.
#' @param fit (Optional) Fitted `locoMDS` object corresponding to the
#'   embeddings in `y`. Only used for warning messages about uniqueness of
#'   points.
#'
#' @returns A data frame with columns:
#' - `k`: Neighborhood size.
#' - `lcmc`: Local Continuity Meta-Criterion value for the given neighborhood
#'   size.
#' - `adjusted_lcmc`: Adjusted LCMC value, accounting for the neighborhood size.
#'
#' @examples
#' # generate example data
#' n <- 100
#' p1 <- 4
#' p2 <- 2
#' x <- matrix(rnorm(n * p1), nrow = n, ncol = p1)
#' y <- matrix(rnorm(n * p2), nrow = n, ncol = p2)
#'
#' # evaluate LCMC for neighborhood sizes 5, 10, and 20
#' lcmc_out <- lcmc(x, y, ks = c(5, 10, 20))
#'
#' @export
lcmc <- function(x, y, ks, fit = NULL) {
  n_samples <- nrow(x)
  k_max <- max(ks)
  y <- check_consensus_embeddings(y, fit = fit)
  x_knns <- FNN::get.knn(x, k = k_max)$nn.index
  y_knns <- FNN::get.knn(y, k = k_max)$nn.index
  result <- purrr::map(
    ks,
    function(k) {
      lcmc_val <- purrr::map_dbl(
        1:n_samples,
        ~ sum(x_knns[.x, 1:k] %in% y_knns[.x, 1:k]) / k
      ) |>
        mean()
      data.frame(
        k = k,
        lcmc = lcmc_val,
        adjusted_lcmc = lcmc_val - (k / (n_samples - 1))
      )
    }
  ) |>
    dplyr::bind_rows()
  return(result)
}


#' Check consensus embeddings for uniqueness of points
#'
#' @inheritParams lcmc
#'
#' @returns Original consensus embeddings if sufficiently unique; otherwise,
#'   the rounded consensus embeddings with a warning.
#' @keywords internal
check_consensus_embeddings <- function(y, fit = NULL) {
  y_rounded <- scale(y, center = TRUE, scale = TRUE) |>
    round(digits = 3)
  y_nunique_per_col <- apply(y_rounded, 2, function(col) length(unique(col)))
  if (any(y_nunique_per_col < 5)) {
    y <- y_rounded
    if (!is.null(fit)) {
      warning(
        paste0(
          "The new consensus space has very few (", min(y_nunique_per_col),
          ") unique points. LCMC results may be unreliable for tau = ",
          fit$tau, " and percentile = ", fit$percentile, "."
        )
      )
    } else {
      warning(
        paste0(
          "The new consensus space has very few (", min(y_nunique_per_col),
          ") unique points. LCMC results may be unreliable."
        )
      )
    }
  }
  return(y)
}
