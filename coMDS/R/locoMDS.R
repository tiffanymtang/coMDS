#' Local Consensus Multidimensional Scaling (LoCoMDS)
#'
#' @inheritParams coMDS
#' @param tau Regularization parameter. Can be a scalar or a vector.
#' @param percentile Percentile used to define local neighborhoods. Larger
#'   percentiles result in larger neighborhoods. Can be a scalar or a vector.
#'
#' @returns A list containing the following components:
#' - `gspace`: Consensus embedding matrix of size `n x ndim`, where `n` is the
#'   number of samples.
#' - `cweights`: A list of learned weight matrices for each input embedding/
#'   distance matrix.
#' - `stress`: The final stress value of the CoMDS solution.
#' - `rss`: Residual sum of squares.
#' - `spp`: Stress per point (in percent).
#' - `spps`: Stress per point per input embedding/distance matrix (in percent,
#'   conditional on the input embedding/distance matrix).
#' - `sps`: Stress per input embedding/distance matrix (in percent).
#' - `sps_norm`: Stress per input embedding/distance matrix (in percent),
#'   normalized by sum of source weights.
#' - `ndim`: Number of dimensions in consensus embedding.
#' - `model`: Type of smacof model.
#' - `tau`: The regularization parameter used in the LoCoMDS algorithm.
#' - `percentile`: The percentile used to define local neighborhoods.
#' - `niter`: Number of iterations run.
#' - `nobj`: Number of points/samples.
#' - `type`: MDS type used (e.g., "interval", "ratio", "ordinal", or "mspline").
#' - `constraint`: Type of weight constraint used in the CoMDS algorithm.
#' - `call`: The matched call to the function.
#' If multiple hyperparameter combinations are provided, the output will be a
#' list of such LoCoMDS objects, each corresponding to a different
#' combination of `tau` and `percentile`.
#'
#' @examples
#' data(iris)
#' # remove duplicates so that tSNE can run
#' X <- dplyr::distinct(iris[, 1:4])
#'
#' # fit various dimension reduction methods
#' pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
#' tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
#' umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
#' dr_list <- list(
#'   pca = pca_scores,
#'   tsne = tsne_scores,
#'   umap = umap_scores
#' )
#'
#' # fit LoCoMDS using dimension reduction embeddings directly as input
#' locomds_out <- locoMDS(
#'   embed_list = dr_list, ndim = 2, tau = 0.1, percentile = 0.5
#' )
#'
#' # fit LoCoMDS using distance matrices as input
#' dist_list <- purrr::map(dr_list, ~ dist(.x, method = "euclidean"))
#' locomds_out <- locoMDS(
#'   D_list = dist_list, ndim = 2, tau = 0.1, percentile = 0.5
#' )
#'
#' # can re-run with larger number of iterations or smaller tolerance if desired
#' locomds_out <- locoMDS(
#'   D_list = dist_list, ndim = 2, tau = 0.1, percentile = 0.5,
#'   itmax = 500, eps = 1e-8
#' )
#'
#' # fit LoCoMDS with multiple possible hyperparameters
#' locomds_out <- locoMDS(
#'   D_list = dist_list, ndim = 2, tau = c(0.01, 0.1), percentile = c(0.5, 0.8)
#' )
#'
#' @export
locoMDS <- function(D_list, embed_list = NULL,
                    tau, percentile, ndim = 2,
                    constraint = c("indscal", "idioscal", "identity"),
                    ...) {
  # setup
  constraint <- match.arg(constraint)
  if (!is.null(embed_list)) {
    D_list <- purrr::map(embed_list, ~ dist(.x, method = "euclidean"))
  }

  # normalize distances for numerical stability
  D_list <- normalize_distances(D_list)

  if ((length(tau) == 1) && (length(percentile) == 1)) {
    # run LoCoMDS for given hyperparameters (tau and percentile)
    result <- localSmacofIndDiff(
      D_list, tau = tau, percentile = percentile, ndim = ndim,
      constraint = constraint, ...
    )
  } else {
    # run LoCoMDS for each combination of hyperparameters (tau and percentile)
    hyperparam_grid <- expand.grid(tau = tau, percentile = percentile)
    result <- furrr::future_map2(
      hyperparam_grid$tau, hyperparam_grid$percentile,
      function(tau, percentile) {
        localSmacofIndDiff(
          D_list, tau = tau, percentile = percentile, ndim = ndim,
          constraint = constraint, ...
        )
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  }
  return(result)
}


#' Tune hyperparameters for LoCoMDS
#'
#' @param obj Output of [locoMDS()].
#' @param data Original data used to evaluate LCMC.
#' @param ks Vector of neighborhood sizes to evaluate LCMC. If `NULL`, defaults to
#'   a range of sizes between 1 and 70% of the number of samples in `data`.
#' @param make_plots Logical indicating whether to create a plot of LCMC values.
#'
#' @returns A list with two components:
#' - `lcmc`: A data frame containing the LCMC values for each neighborhood size
#'   and hyperparameter combination.
#' - `plots`: A ggplot object showing the LCMC values across neighborhood sizes
#'   and hyperparameters.
#'
#' @examples
#' data(iris)
#' # remove duplicates so that tSNE can run
#' X <- dplyr::distinct(iris[, 1:4])
#'
#' # fit various dimension reduction methods
#' pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
#' tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
#' umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
#' dr_list <- list(
#'   pca = pca_scores,
#'   tsne = tsne_scores,
#'   umap = umap_scores
#' )
#'
#' # fit LoCoMDS with multiple possible hyperparameters
#' locomds_fits_all <- locoMDS(
#'   embed_list = dr_list, ndim = 2, tau = c(0.01, 0.1), percentile = c(0.5, 0.8)
#' )
#'
#' # tune LoCoMDS hyperparameters
#' locomds_tuning <- tune_locoMDS(locomds_fits_all, data = X)
#'
#' @export
tune_locoMDS <- function(obj, data, ks = NULL, make_plots = TRUE) {
  if (inherits(obj, "locoMDS")) {
    warning(
      "locoMDS() was run with a single set of hyperparameters. ",
      "No tuning is necessary. ",
      "Returning the original locoMDS object."
    )
    return(obj)
  } else if (!is.list(obj)) {
    stop("Input `obj` should be the output of `locoMDS()`.")
  }
  if (is.null(ks)) {
    ks <- c(1, 2, 5, round(10^seq(1, log10(nrow(data) * 0.7), length.out = 7)))
  }
  # compute lcmc
  lcmc_result <- purrr::map(
    obj,
    ~ lcmc(x = data, y = .x$gspace, ks = ks, fit = .x) |>
      dplyr::mutate(
        tau = .x$tau,
        percentile = .x$percentile
      )
  ) |>
    dplyr::bind_rows() |>
    dplyr::select(k, tau, percentile, lcmc, adjusted_lcmc)
  # make lcmc plots
  plt <- NULL
  if (make_plots) {
    # get best percentile for each tau
    lcmc_result_taus <- lcmc_result |>
      dplyr::group_by(tau, k) |>
      dplyr::slice_max(adjusted_lcmc, n = 1) |>
      dplyr::group_by(tau, percentile) |>
      dplyr::summarise(
        freq = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::group_by(tau) |>
      dplyr::slice_max(freq, n = 1) |>
      dplyr::select(-freq)

    # get best tau for each percentile
    lcmc_result_percentiles <- lcmc_result |>
      dplyr::group_by(percentile, k) |>
      dplyr::slice_max(adjusted_lcmc, n = 1) |>
      dplyr::group_by(percentile, tau) |>
      dplyr::summarise(
        freq = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::group_by(percentile) |>
      dplyr::slice_max(freq, n = 1) |>
      dplyr::select(-freq)

    # make plots:
    # (i) full plot with all taus and percentiles
    # (ii) abridged plot with only the best percentile for each tau
    # (iii) abridged plot with only the best tau for each percentile
    plt_df <- dplyr::bind_rows(
      lcmc_result |>
        dplyr::mutate(.mode = "All"),
      lcmc_result |>
        dplyr::inner_join(lcmc_result_taus, by = c("tau", "percentile")) |>
        dplyr::mutate(.mode = "Taus with best percentile"),
      lcmc_result |>
        dplyr::inner_join(lcmc_result_percentiles, by = c("tau", "percentile")) |>
        dplyr::mutate(.mode = "Percentiles with best tau")
    )
    plt <- plot_locoMDS_lcmc(plt_df) +
      ggplot2::facet_grid(~ .mode)
  }
  return(list(
    lcmc = lcmc_result,
    plot = plt
  ))
}


#' Extract LoCoMDS fit for particular tau and percentile
#'
#' @inheritParams locoMDS
#' @inheritParams tune_locoMDS
#'
#' @returns The locoMDS object corresponding to the provided tau and percentile
#'
#' @examples
#' data(iris)
#' # remove duplicates so that tSNE can run
#' X <- dplyr::distinct(iris[, 1:4])
#'
#' # fit various dimension reduction methods
#' pca_scores <- prcomp(X, center = TRUE, scale = TRUE)$x
#' tsne_scores <- Rtsne::Rtsne(X, dims = 2, perplexity = 30, verbose = FALSE)$Y
#' umap_scores <- umap::umap(X, n_components = 2, verbose = FALSE)$layout
#' dr_list <- list(
#'   pca = pca_scores,
#'   tsne = tsne_scores,
#'   umap = umap_scores
#' )
#'
#' # fit LoCoMDS with multiple possible hyperparameters
#' locomds_fits_all <- locoMDS(
#'   embed_list = dr_list, ndim = 2, tau = c(0.01, 0.1), percentile = c(0.5, 0.8)
#' )
#'
#' # extract LoCoMDS with specific hyperparameters
#' locomds_fit <- extract_locoMDS(locomds_fits_all, tau = 0.1, percentile = 0.5)
#'
#' @export
extract_locoMDS <- function(obj, tau, percentile) {
  if (inherits(obj, "locoMDS")) {
    warning(
      "locoMDS() was run with a single set of hyperparameters. ",
      "Returning the original locoMDS object."
    )
    return(obj)
  } else if (!is.list(obj)) {
    stop("Input `obj` should be the output of `locoMDS()`.")
  }
  tuned_idx <- which(
    sapply(
      obj,
      function(x) {
        isTRUE(all.equal(x$tau, tau)) && isTRUE(all.equal(x$percentile, percentile))
      }
    )
  )
  return(obj[[tuned_idx]])
}
