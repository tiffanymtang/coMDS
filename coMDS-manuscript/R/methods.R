#' Fit individual dimension reduction methods
#'
#' @param data A matrix or data frame containing covariates data
#' @param methods_list Named list of dimension reduction method functions to fit
fit_dimension_reduction_methods <- function(data, methods_list) {
  dist_data <- dist(data)
  embeddings_out <- purrr::map(
    methods_list,
    function(f) {
      embedding_out <- f(data = data, dist_data = dist_data)
      rownames(embedding_out) <- 1:nrow(embedding_out)
      colnames(embedding_out) <- paste("Component", 1:ncol(embedding_out))
      return(embedding_out)
    }
  )
  return(embeddings_out)
}


#' Fit Consensus Multidimensional Scaling (CoMDS)
#'
#' @param dist_list A list of distance matrices.
#' @param data_list A list of dimension reduction embeddings. Only used (and
#'   required) if dist_list is not provided.
#' @param ndim Number of dimension reduction components
#' @param ... Additional arguments to be passed to coMDS::coMDS
#'
#' @returns A list of two elements:
#' - embeddings: a data frame of the CoMDS embedding
#' - fit: the CoMDS model fit
fit_comds <- function(dist_list, data_list = NULL, ndim = 2, ...) {
  if (missing(dist_list) && is.null(data_list)) {
    stop("Either 'dist_list' or 'data_list' must be provided.")
  } else if (missing(dist_list)) {
    dist_list <- purrr::map(data_list, ~ dist(.x))
  }
  comds_out <- coMDS::coMDS(D_list = dist_list, ndim = ndim, ...)
  embedding_out <- as.data.frame(comds_out$gspace) |>
    setNames(paste("Component", 1:ndim))
  return(list(
    embeddings = embedding_out,
    fit = comds_out
  ))
}


#' Fit Local Consensus Multidimensional Scaling (LoCoMDS)
#'
#' @param dist_list A list of distance matrices.
#' @param data_list A list of dimension reduction embeddings. Only used (and
#'   required) if dist_list is not provided.
#' @param data The original data matrix; for tuning LoCoMDS hyperparameters.
#' @param ndim Number of dimension reduction components.
#' @param taus A vector of tau hyperparameter values to try.
#' @param percentiles A vector of percentile hyperparameter values to try.
#' @param ... Additional arguments to be passed to coMDS::locoMDS
#'
#' @returns A list of two elements:
#' - embeddings: a data frame of the CoMDS embedding
#' - fit: the LoCoMDS model fit
#' - tuning: the LoCoMDS hyperparameter tuning results
#' - best_fit: the LoCoMDS model fit with the best hyperparameters
fit_locomds <- function(dist_list, data_list = NULL, data, ndim = 2,
                        taus = c(100, 20, 10, 1, 0.1, 0.05, 0.01, 0.005, 0.001),
                        percentiles = seq(0.1, 0.9, 0.2), min_distinct_pts = NULL, ...) {
  if (missing(dist_list) && is.null(data_list)) {
    stop("Either 'dist_list' or 'data_list' must be provided.")
  } else if (missing(dist_list)) {
    dist_list <- purrr::map(data_list, ~ dist(.x))
  }
  locomds_out <- coMDS::locoMDS(
    D_list = dist_list, ndim = ndim, tau = taus, percentile = percentiles, ...
  )

  locomds_tuned_out <- coMDS::tune_locoMDS(locomds_out, data)

  best_params <- locomds_tuned_out$lcmc |>
  dplyr::group_by(tau, percentile) |>
  dplyr::slice_max(adjusted_lcmc, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  (\(df) {
    k_prime <- df |>
      dplyr::count(k) |>
      dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
      dplyr::pull(k)

    locomds_tuned_out$lcmc |>
      dplyr::filter(k < k_prime) |>
      dplyr::group_by(tau, percentile) |>
      dplyr::summarise(
        max_adjusted_lcmc = max(adjusted_lcmc),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(max_adjusted_lcmc))  # keep original order among ties
  })()

  if (is.null(min_distinct_pts)){
      min_distinct_pts <- dist_list[[1]] |> nrow() * 0.01
      print(min_distinct_pts)
  }

  for (i in seq_len(nrow(best_params))) {
    this_emb <- coMDS::extract_locoMDS(
      locomds_out,
      tau        = best_params$tau[i],
      percentile = best_params$percentile[i]
    )

    this_emb_rounded <- this_emb$gspace |>
      scale(center = TRUE, scale = TRUE) |>
      round(3)

    nunique_per_col <- apply(this_emb_rounded, 2, function(x) length(unique(x)))

    if (all(nunique_per_col >= min_distinct_pts)) {
      stable_best <- best_params[i, ]
      break
    }
  }

  best_locomds_out <- coMDS::extract_locoMDS(
    locomds_out, stable_best$tau[1], stable_best$percentile[1]
  )
  embedding_out <- as.data.frame(best_locomds_out$gspace) |>
    setNames(paste("Component", 1:ndim))
  return(list(
    embeddings = embedding_out,
    fit = locomds_out,
    tuning = locomds_tuned_out,
    best_fit = best_locomds_out
  ))
}


#' Fit Meta-Spec
#'
#' @param data_list A list of dimension reduction embeddings
#' @param methods Methods for meta-visualization; one or more of "kPCA" and
#'   "UMAP"
#' @param ndim Number of dimension reduction components for meta-visualization
#'
#' @returns A list of two elements:
#' - embeddings: a list of meta-visualization embeddings
#' - fit: the metaspec model fit
fit_metaspec <- function(data_list, methods = c("kPCA", "UMAP"), ndim = 2) {
  methods <- match.arg(methods, several.ok = TRUE)
  ensemble_out <- ensemble_metaspec(data_list)
  meta_dist_mat <- ensemble_out$ensemble_dist_mat

  embeddings_out <- list()
  if ("kPCA" %in% methods) {
    # Converts the distance matrix into a Gaussian kernel matrix, 
    # then Uses the median distance as the kernel bandwidth
    K_mat <- exp(-as.matrix(meta_dist_mat)^2 / quantile(meta_dist_mat, 0.5)^2)
    # Compute eigenvalues and eigenvectors of the kernel matrix
    eigen_K_mat <- eigen(K_mat)
    # skip the first eigenvector (often corresponds to trivial/invariant direction in kernel PCA)
    # scale each eigenvector by the square root (or here: by the eigenvalue directly)
    embeddings_out[["kPCA"]] <- eigen_K_mat$vectors[, 2:(ndim + 1)] %*%
      diag(eigen_K_mat$values[2:(ndim + 1)])
  }
  if ("UMAP" %in% methods) {
    embeddings_out[["UMAP"]] <- uwot::umap(
      as.dist(meta_dist_mat),  n_neighbors = 30
    )
  }
  embeddings_out <- purrr::map(
    embeddings_out,
    function(x) {
      as.data.frame(x) |> setNames(paste("Component", 1:ndim))
    }
  )
  return(list(
    embeddings = embeddings_out,
    fit = ensemble_out
  ))
}


#' Fit Multi-SNE
#'
#' @param data_list A list of dimension reduction embeddings.
#' @param perplexity_vec A vector of perplexity values.
#' @param ndim Number of dimension reduction components.
#' @param ... Additional arguments to be passed to multiSNE::multiSNE
#'
#' @returns A list of two elements:
#' - embeddings: a list of multi-SNE embeddings
#' - fit: a list of the multi-SNE model fit for each perplexity
fit_multisne <- function(data_list, perplexity_vec = 30, ndim = 2, ...) {
  fit_list <- purrr::map(
    perplexity_vec,
    function(perplexity) {
      multiSNE::multiSNE(X = data_list, k = ndim, perplexity = perplexity, ...)
    }
  ) |> setNames(paste0("perplexity_", perplexity_vec)) 

  embeddings_out <- purrr::map(
    fit_list,
    function(fit) {
      as.data.frame(fit$Y) |> setNames(paste("Component", 1:ndim))
    }
  )
  return(list(
    embeddings = embeddings_out,
    fit = fit_list  ))
}
