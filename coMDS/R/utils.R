#' Normalize list of distance matrices
#'
#' @description
#' This function normalizes a list of distance matrices by their top singular
#' value so that after normalization, each distance matrix has a top singular
#' value of 1. This ensures that the distance matrices are on similar scales,
#' making the CoMDS optimization problem more stable and easier to optimize.
#'
#' @param D_list List of distance matrices.
#'
#' @returns A list of normalized distance matrices (each of class `dist`).
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
#' # convert embeddings to distance matrices
#' dist_list <- purrr::map(dr_list, ~ dist(.x, method = "euclidean"))
#' # normalize distance matrices
#' dist_list_norm <- normalize_distances(dist_list)
#'
#' @export
normalize_distances <- function(D_list) {
  D_list <- purrr::map(
    D_list,
    function(D) {
      D_mat <- as.matrix(D)
      if (any(is.na(D_mat))) {
        diag(D_mat) <- NA
        na_rows <- which(apply(D_mat, 1, function(x) all(is.na(x))))
        D_mat <- D_mat[-na_rows, -na_rows]
        diag(D_mat) <- 0
      }
      D_svd <- svd(D_mat, nu = 1, nv = 1)
      D_norm <- D / D_svd$d[1]
      D_norm[D_norm < 0] <- 0
      return(as.dist(D_norm))
    }
  )
  return(D_list)
}
