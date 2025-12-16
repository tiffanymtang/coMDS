#' Fit Principal Component Analysis (PCA)
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
fit_pca <- function(data, ndim = 2, ...) {
  pca_out <- rARPACK::svds(as.matrix(data), k = ndim)
  return(pca_out$u)
}

#' Fit classical Multidimensional Scaling (MDS)
#'
#' @param dist_data A distance matrix
#' @param ndim Number of dimension reduction components
fit_mds <- function(dist_data, ndim = 2, ...) {
  mds_out <- cmdscale(dist_data, k = ndim)
  return(mds_out)
}

#' Fit non-metric Multidimensional Scaling (Non-metric MDS)
#'
#' @param dist_data A distance matrix
#' @param ndim Number of dimension reduction components
fit_imds <- function(dist_data, ndim = 2, verbose = FALSE, ...) {
  imds_out <- MASS::isoMDS(dist_data, k = ndim, trace = verbose)
  return(imds_out$points)
}

#' Fit Sammon's nonlinear mapping
#'
#' @param dist_data A distance matrix
#' @param ndim Number of dimension reduction components
fit_sammon <- function(dist_data, ndim = 2, verbose = FALSE, ...) {
  sammon_out <- MASS::sammon(dist_data, k = ndim, trace = verbose)
  return(sammon_out$points)
}

#' Fit Locally Linear Embedding (LLE)
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
#' @param n_neighbors Number of neighbors
#' @param reg Regularization method; one of 1, 2, or 3. See [lle::lle()] for
#'   details.
fit_lle <- function(data, ndim = 2, n_neighbors = 20, reg = 2, verbose = FALSE,
                    ...) {
  if (verbose) {
    lle_out <- lle::lle(as.matrix(data), m = ndim, k = n_neighbors, reg = reg)
  } else {
    sink(tempfile())
    on.exit(sink())
    lle_out <- lle::lle(as.matrix(data), m = ndim, k = n_neighbors, reg = reg)
  }
  return(lle_out$Y)
}

#' Fit Hessian Locally Linear Embedding (HLLE)
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
#' @param n_neighbors Number of neighbors
fit_hlle <- function(data, ndim = 2, n_neighbors = 20, verbose = FALSE, ...) {
  if (verbose) {
    .mute <- character(0)
  } else {
    .mute <- c("message", "output")
  }
  hlle_out <- dimRed::embed(
    data, "HLLE", ndim = ndim, knn = n_neighbors, .mute = .mute
  )
  return(hlle_out@data@data)
}

#' Fit Isomap
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
#' @param n_neighbors Number of neighbors
fit_isomap <- function(data, ndim = 2, n_neighbors = 20, verbose = FALSE, ...) {
  if (verbose) {
    .mute <- character(0)
  } else {
    .mute <- c("message", "output")
  }
  isomap_out <- dimRed::embed(
    data, "Isomap", ndim = ndim, knn = n_neighbors, .mute = .mute
  )
  return(isomap_out@data@data)
}

#' Fit kernel PCA (kPCA)
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
#' @param sigma Bandwidth parameter for the RBF kernel
fit_kpca <- function(data, ndim = 2, sigma = 0.01, verbose = FALSE, ...) {
  if (verbose) {
    .mute <- character(0)
  } else {
    .mute <- c("message", "output")
  }
  kpca_out <- dimRed::embed(
    data, "kPCA", ndim = ndim, kpar = list(sigma = sigma), .mute = .mute
  )
  return(kpca_out@data@data)
}

#' Fit Laplacian Eigenmap (LEIM)
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
fit_leim <- function(data, ndim = 2, verbose = FALSE, ...) {
  if (verbose) {
    .mute <- character(0)
  } else {
    .mute <- c("message", "output")
  }
  leim_out <- dimRed::embed(
    data, "LaplacianEigenmaps", ndim = ndim, .mute = .mute
  )
  return(leim_out@data@data)
}

#' Fit UMAP
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
#' @param n_neighbors Number of neighbors
#' @param min_dist Minimum distance between points in the low-dimensional space
fit_umap <- function(data, ndim = 2, n_neighbors = 30, min_dist = 0.01, ...) {
  umap_out <- uwot::umap(
    data, n_components = ndim, n_neighbors = n_neighbors, min_dist = min_dist
  )
  return(umap_out)
}

#' Fit tSNE
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
#' @param perplexity Perplexity parameter for tSNE
fit_tsne <- function(data, ndim = 2, perplexity = 30, verbose = FALSE, ...) {
  if (verbose) {
    .mute <- character(0)
  } else {
    .mute <- c("message", "output")
  }
  tsne_out <- dimRed::embed(
    data, "tSNE", ndim = ndim, perplexity = perplexity, .mute = .mute
  )
  return(tsne_out@data@data)
}

#' Fit PHATE
#'
#' @param data A matrix or data frame containing covariates data
#' @param ndim Number of dimension reduction components
#' @param n_neighbors Number of neighbors
fit_phate <- function(data, ndim = 2, n_neighbors = 30, verbose = FALSE, ...) {
  phate_out <- phateR::phate(
    as.matrix(data), ndim = ndim, knn = n_neighbors, verbose = verbose, seed = 331
  )
  return(phate_out$embedding)
}
